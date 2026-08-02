"""Gene-drug pair references.

We store identifiers and URLs only, never guideline prose. Two reasons:
guidelines are revised, so a stored copy is a stale copy; and redistributing
CPIC/PharmGKB content commercially may require permission we have not
obtained. Linking sidesteps both.

The table is validated on load rather than trusted. A pair silently dropped
because of a typo in the JSON does not surface as a loading error -- it
surfaces as `no_guidance_for_pair` for a drug CPIC does publish for, which is
the exact false negative this project exists to avoid. So every entry must be
complete, non-blank, https, and unique, or nothing loads at all.
"""

from __future__ import annotations

import json
import unicodedata
from dataclasses import dataclass, fields
from pathlib import Path
from urllib.parse import urlsplit


class GuidelineTableError(ValueError):
    """The gene-drug pair table cannot be trusted, so it is not loaded.

    Raised instead of skipping the bad entries. A partially loaded table
    answers "CPIC publishes nothing for this drug" for pairs it merely failed
    to parse, and that answer is indistinguishable from the truthful one.
    """


@dataclass(frozen=True)
class GuidelineRef:
    """A pointer to published guidance for one gene-drug pair."""

    gene: str
    drug: str
    cpic_pair_id: str
    url: str


_FIELD_NAMES = tuple(f.name for f in fields(GuidelineRef))

# Identifiers and links are short. A field this long is prose, and prose does
# not belong in this table (see the module docstring).
_MAX_FIELD_LENGTH = 200

# The only host we will cite. The entire output of this tool is a citation, so
# a row pointing anywhere else is a confident answer sourced from somewhere we
# never vetted -- indistinguishable, to a reader, from a real one.
_CITATION_HOST = "cpicpgx.org"


def normalize_drug(drug: str) -> str:
    """The comparison form of a drug name: NFKC-normalized, stripped, casefolded.

    Spelled once and used by every lookup so that a query and the table can
    never be normalized differently. Rejects a blank name rather than
    returning "" -- a blank query that comes back "no guidance for this drug"
    is a fabricated negative about a drug nobody named.

    NFKC folds the compatibility forms a real paste produces -- fullwidth
    "ＷＡＲＦＡＲＩＮ" from an IME or a PDF, a non-breaking space -- onto their
    ASCII equivalents. Without it those queries fall through to
    `no_guidance_for_pair`, which is the safe direction but still a wrong
    answer nobody would question. NFKC does *not* fold homoglyphs from other
    scripts (Cyrillic "а" stays distinct from Latin "a"), so a spoofed name
    still fails closed, which is correct: we would rather refuse to recognize
    a name than match the wrong drug.
    """
    if not isinstance(drug, str):
        raise TypeError(f"drug name must be a string, got {type(drug).__name__}")
    # NFKC first: it turns U+00A0 and friends into plain spaces, so strip()
    # only works on the whole of the whitespace after normalizing.
    needle = unicodedata.normalize("NFKC", drug).strip().casefold()
    if not needle:
        raise ValueError(
            f"refusing to look up a blank drug name ({drug!r}): a blank query "
            f"answered 'no guidance' would be a negative finding about a drug "
            f"that was never named"
        )
    return needle


def _entry_to_ref(index: int, entry: object) -> GuidelineRef:
    """Validate one raw table entry and build its GuidelineRef."""
    if not isinstance(entry, dict):
        raise GuidelineTableError(
            f"entry {index} is not an object (got {type(entry).__name__}); the "
            f"pair table must be a list of objects"
        )

    missing = [name for name in _FIELD_NAMES if name not in entry]
    if missing:
        raise GuidelineTableError(
            f"entry {index} is missing required field(s) {', '.join(missing)}"
        )
    # An unexpected key is most likely embedded guideline prose or a renamed
    # field; either way GuidelineRef(**entry) would raise a bare TypeError, so
    # this says what is actually wrong.
    unexpected = sorted(set(entry) - set(_FIELD_NAMES))
    if unexpected:
        raise GuidelineTableError(
            f"entry {index} has unexpected field(s) {', '.join(unexpected)}; this "
            f"table carries identifiers and links only, never guideline text"
        )

    for name in _FIELD_NAMES:
        value = entry[name]
        if not isinstance(value, str):
            raise GuidelineTableError(
                f"entry {index} field {name!r} must be a string, got "
                f"{type(value).__name__}"
            )
        if not value.strip():
            raise GuidelineTableError(f"entry {index} field {name!r} is blank")
        if len(value) > _MAX_FIELD_LENGTH:
            raise GuidelineTableError(
                f"entry {index} field {name!r} is {len(value)} characters; this "
                f"table stores identifiers and links, not guideline prose"
            )

    url = entry["url"]
    if not url.startswith("https://"):
        raise GuidelineTableError(
            f"entry {index} url {url!r} is not https; guidance must be cited over "
            f"a channel that cannot be rewritten in transit"
        )

    # A row whose citation does not match its own gene is the worst failure this
    # table has, because it is silent: the query answers guidance_found and
    # cites a real-looking guideline for a different gene. Nothing downstream
    # can catch it -- the citation IS the answer -- so it has to be caught here.
    gene = entry["gene"]
    if gene.strip().casefold() not in entry["cpic_pair_id"].casefold():
        raise GuidelineTableError(
            f"entry {index} cites cpic_pair_id {entry['cpic_pair_id']!r}, which "
            f"does not name its own gene {gene!r}; a pair id belonging to another "
            f"gene would be reported as guidance for this one"
        )

    # Only the gene is checked against the pair id, never the drug: CPIC names
    # some guidelines after a drug class rather than a member, so DPYD's real
    # pair id could legitimately be "DPYD-fluoropyrimidines" for a row whose
    # drug is "capecitabine". Requiring the drug to appear would reject correct
    # future rows, and a rejected row loads as nothing at all.
    host = urlsplit(url).hostname
    if host is None or (
        host != _CITATION_HOST and not host.endswith(f".{_CITATION_HOST}")
    ):
        raise GuidelineTableError(
            f"entry {index} url {url!r} is not on {_CITATION_HOST} (host "
            f"{host!r}); this table cites CPIC and nothing else, so a link "
            f"elsewhere would be presented as CPIC guidance"
        )

    return GuidelineRef(**entry)


def load_pairs(path: Path) -> list[GuidelineRef]:
    """Load the gene-drug pair reference table, or raise.

    Never returns a partial table and never returns an empty one: an empty
    table makes every drug look like a drug CPIC does not cover.
    """
    path = Path(path)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as err:
        raise GuidelineTableError(
            f"could not read the gene-drug pair table {path}: {err}"
        ) from err

    if not isinstance(raw, list):
        raise GuidelineTableError(
            f"gene-drug pair table {path} must be a JSON list, got "
            f"{type(raw).__name__}"
        )
    if not raw:
        raise GuidelineTableError(
            f"gene-drug pair table {path} contains no gene-drug pairs; an empty "
            f"table would report every drug as one CPIC does not publish for"
        )

    pairs = [_entry_to_ref(index, entry) for index, entry in enumerate(raw)]

    seen: dict[tuple[str, str], int] = {}
    seen_ids: dict[str, int] = {}
    for index, pair in enumerate(pairs):
        key = (pair.gene, normalize_drug(pair.drug))
        if key in seen:
            raise GuidelineTableError(
                f"entry {index} is a duplicate gene-drug pair "
                f"{pair.gene}/{pair.drug} (first seen at entry {seen[key]}); a "
                f"duplicated pair would be reported twice for one query"
            )
        seen[key] = index
        if pair.cpic_pair_id in seen_ids:
            raise GuidelineTableError(
                f"entry {index} reuses duplicate cpic_pair_id "
                f"{pair.cpic_pair_id!r} (first seen at entry "
                f"{seen_ids[pair.cpic_pair_id]}); Task 8 keys guideline changes "
                f"on this id, so it must identify one pair"
            )
        seen_ids[pair.cpic_pair_id] = index

    return pairs


def find_pairs_for_drug(drug: str, pairs: list[GuidelineRef]) -> list[GuidelineRef]:
    """Every gene-drug pair matching a drug name, case-insensitively.

    Matching ignores case and surrounding whitespace on both sides. A user who
    types "Clopidogrel" and a user who types "clopidogrel" must not get
    different answers; a lookup that misses is silently indistinguishable from
    a drug CPIC does not cover.
    """
    needle = normalize_drug(drug)
    return [p for p in pairs if normalize_drug(p.drug) == needle]
