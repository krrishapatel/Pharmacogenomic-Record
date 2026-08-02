"""Answer drug queries against a stored record.

The load-bearing rule of this module: absence of guidance and absence of data
are different answers and must never collapse into one. A gene a consumer
array never covered carries NO information, and presenting that as "no
interaction found" is the most dangerous thing this system could do.

Hence three outcomes, always distinguishable:

  guidance_found       gene called, CPIC publishes guidance for this pair
  no_guidance_for_pair gene called, CPIC publishes nothing for this drug
  cannot_assess        gene not covered, indeterminate, or absent -- unknown

Two further rules follow from that one:

* What comes back is a *reference*, never a recommendation. Every explanation
  names the pair and cites the URL; none of them states a dose, a drug choice,
  or what anyone should do. The guideline text lives at CPIC, where it is kept
  current, and is not copied here.
* Silence is the failure mode with consequences, so nothing here answers with
  an empty list, and nothing swallows an error. A `CorruptRecordError` from the
  store propagates: a database integrity failure must not be reported as a mere
  coverage gap, because "we cannot read this record" and "the array missed this
  gene" call for completely different actions.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from pgxrecord.caller import CALLED, GeneCall
from pgxrecord.guidelines import GuidelineRef, find_pairs_for_drug, normalize_drug
from pgxrecord.store import CorruptRecordError, RecordStore

GUIDANCE_FOUND = "guidance_found"
NO_GUIDANCE_FOR_PAIR = "no_guidance_for_pair"
CANNOT_ASSESS = "cannot_assess"

OUTCOMES = (GUIDANCE_FOUND, NO_GUIDANCE_FOR_PAIR, CANNOT_ASSESS)

# How reassuring each outcome is, least first. Used to summarize a multi-gene
# answer: the summary may never be more confident than its weakest component.
# cannot_assess outranks both because an unknown gene invalidates any claim of
# completeness; guidance_found outranks no_guidance_for_pair because a real
# pair must not be hidden behind a gene CPIC happens not to publish for.
_SEVERITY = {
    CANNOT_ASSESS: 2,
    GUIDANCE_FOUND: 1,
    NO_GUIDANCE_FOR_PAIR: 0,
}

# Repeated in every cannot_assess explanation. The wording is deliberate: it
# names the state as missing data, so no consumer -- human or string-matching
# -- can read the answer as a negative finding.
_NOT_A_NEGATIVE = "That is missing data, not absence of an interaction."


@dataclass(frozen=True)
class QueryResult:
    """One gene's answer for a drug query."""

    outcome: str
    gene: str | None
    phenotype: str | None
    guideline: GuidelineRef | None
    explanation: str


def _calls_by_gene(calls: Iterable[GeneCall]) -> dict[str, GeneCall]:
    """Index gene calls by gene symbol, refusing to resolve a conflict.

    A plain dict comprehension is last-write-wins, which would let a second
    `called` row for a gene mask a `not_covered` one -- the collapse this whole
    module exists to prevent, arriving through the back door. The store's
    (record_id, gene) primary key makes this unreachable via `append()`, so a
    duplicate means the file was written by something else, and the honest
    answer is to refuse rather than pick.
    """
    by_gene: dict[str, GeneCall] = {}
    for call in calls:
        previous = by_gene.get(call.gene)
        if previous is not None:
            raise CorruptRecordError(
                f"record stores two calls for {call.gene} "
                f"(coverage {previous.coverage!r} and {call.coverage!r}); "
                f"refusing to choose between them, because picking the called "
                f"one would hide a gene that was never covered"
            )
        by_gene[call.gene] = call
    return by_gene


def _no_pair_result(drug: str) -> QueryResult:
    """The explicit negative for a drug this table has no pair for."""
    return QueryResult(
        outcome=NO_GUIDANCE_FOR_PAIR,
        gene=None,
        phenotype=None,
        guideline=None,
        explanation=(
            f"This reference table lists no CPIC gene-drug pair for {drug!r}. "
            f"The table is a curated subset of CPIC's pairs, so this states only "
            f"that {drug!r} is absent from it."
        ),
    )


def _unassessable(pair: GuidelineRef, reason: str) -> QueryResult:
    """A cannot_assess result. Phenotype is always None; there is not one."""
    return QueryResult(
        outcome=CANNOT_ASSESS,
        gene=pair.gene,
        phenotype=None,
        guideline=pair,
        explanation=(
            f"Cannot assess {pair.gene} for pair {pair.cpic_pair_id}: {reason} "
            f"{_NOT_A_NEGATIVE} CPIC guidance for this pair is at {pair.url}."
        ),
    )


def _assessed(pair: GuidelineRef, call: GeneCall) -> QueryResult:
    """A guidance_found result: a citation, never a recommendation.

    The phenotype is reported as PharmCAT assigned it, and may legitimately be
    None -- F2, F5, VKORC1, CFTR, IFNL3 and ABCG2 have no metabolizer phenotype
    at all. Only `coverage` decides whether a gene is assessable, so a null
    phenotype here is a called gene and is reported as one.
    """
    if call.phenotype is None:
        phenotype_text = (
            f"PharmCAT assigns this gene no metabolizer phenotype; the diplotype "
            f"is what the guideline keys on."
        )
    else:
        phenotype_text = f"PharmCAT assigned phenotype {call.phenotype!r}."
    return QueryResult(
        outcome=GUIDANCE_FOUND,
        gene=pair.gene,
        phenotype=call.phenotype,
        guideline=pair,
        explanation=(
            f"{pair.gene} diplotype {call.diplotype}. {phenotype_text} "
            f"CPIC publishes guidance for pair {pair.cpic_pair_id}; the guideline "
            f"itself is at {pair.url} and is not reproduced here."
        ),
    )


def query_drug(
    store: RecordStore,
    subject_id: str,
    drug: str,
    pairs: list[GuidelineRef],
) -> list[QueryResult]:
    """Return one result per gene-drug pair relevant to this drug.

    Never returns an empty list. A drug with no pair in the table gets one
    explicit `no_guidance_for_pair` result; every pair the table does list gets
    its own result, so a drug with several relevant genes cannot report the
    covered ones and quietly drop the rest.

    Raises rather than answering when the query itself is meaningless (a blank
    or non-string drug name) or when the stored record cannot be read
    (`CorruptRecordError`, propagated from the store).
    """
    # Validated before the store is touched: a blank drug name answered
    # "no guidance" would be a negative finding about a drug nobody named.
    normalize_drug(drug)
    # An empty table makes every drug in existence come back
    # no_guidance_for_pair. load_pairs already refuses to return one, so this
    # only catches a caller that built the list by hand -- but the failure is
    # silent and uniform, which is exactly the kind that ships.
    if not pairs:
        raise ValueError(
            "refusing to query against an empty gene-drug pair table: every "
            "drug would be reported as one CPIC does not publish for"
        )

    relevant = find_pairs_for_drug(drug, pairs)
    if not relevant:
        return [_no_pair_result(drug.strip())]

    # Deliberately not wrapped in try/except. latest() returns [] only when the
    # subject has no records at all; a record that exists but cannot be read
    # raises CorruptRecordError, and that must reach the caller as an integrity
    # failure rather than be recast as a coverage gap.
    stored = store.latest(subject_id)
    calls = _calls_by_gene(stored)

    results: list[QueryResult] = []
    for pair in relevant:
        call = calls.get(pair.gene)

        if not stored:
            results.append(
                _unassessable(
                    pair,
                    f"no record is stored for subject {subject_id!r}, so this "
                    f"gene's genotype is unknown.",
                )
            )
            continue

        if call is None:
            results.append(
                _unassessable(
                    pair,
                    f"the latest stored record for subject {subject_id!r} holds "
                    f"no call for this gene, so its genotype is unknown.",
                )
            )
            continue

        if call.coverage != CALLED:
            results.append(
                _unassessable(
                    pair,
                    f"the stored coverage state is {call.coverage!r}, so this "
                    f"gene's genotype is unknown.",
                )
            )
            continue

        results.append(_assessed(pair, call))

    return results


def overall_outcome(results: Sequence[QueryResult]) -> str:
    """Summarize a multi-gene answer as its least reassuring component.

    Exists so that a caller rendering one line per drug cannot accidentally
    report `guidance_found` for warfarin when CYP2C9 was called and VKORC1 was
    never covered. A summary that hides an unknown gene is the collapse this
    module forbids, one level up.

    Refuses to summarize an empty list or an unrecognized outcome: both would
    have to be given a default, and every safe default here is a lie.
    """
    if not results:
        raise ValueError(
            "cannot summarize an empty result list; query_drug never returns "
            "one, so an empty list means a caller dropped results"
        )
    unknown = sorted({r.outcome for r in results} - set(OUTCOMES))
    if unknown:
        raise ValueError(
            f"unrecognized outcome(s) {', '.join(repr(u) for u in unknown)}; the "
            f"vocabulary is exactly {OUTCOMES}"
        )
    return max((r.outcome for r in results), key=lambda outcome: _SEVERITY[outcome])
