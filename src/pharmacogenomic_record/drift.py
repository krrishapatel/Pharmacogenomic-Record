"""Report which stored records a guideline revision touches.

This is the reason the record store exists. A batch tool answers "what does
this genotype mean today". A persistent record can answer "whose stored
results changed meaning when the guidance moved" -- including subjects whose
gene was never covered, since new guidance may be the very reason to finally
get proper testing.

Two rules follow from what this report is used for.

* Coverage is never a filter. `called`, `not_covered` and `indeterminate` all
  land in the report, because the question here is "could this revision concern
  this subject", and a gene nobody has looked at is not a gene a revision
  cannot concern. Filtering on coverage would quietly shrink the report to the
  subjects who happen to have data, which is the population least in need of
  being told.
* An omission here is silent. Every other module in this project can answer
  "cannot assess"; this one answers with a list, and a list that is short by one
  line is indistinguishable from a revision that did not concern that person.
  That is why the id comparison is normalized on both sides and why an empty
  pair table is refused rather than reported as "nobody affected".

What comes back is identifiers only -- subject, gene, pair ids. No guideline
prose, no dose, no imperative: the revised guidance lives at CPIC, and this
module says only where to look, never what to do.
"""

from __future__ import annotations

from dataclasses import dataclass

from pharmacogenomic_record.guidelines import (
    GuidelineRef,
    normalize_gene,
    normalize_pair_id,
)
from pharmacogenomic_record.store import RecordStore


@dataclass(frozen=True)
class AffectedRecord:
    """A stored subject-gene pair touched by a guideline change.

    `changed_pair_ids` is a tuple, not a list, and that is a correctness
    choice rather than a stylistic one: a frozen dataclass holding a list is
    only shallowly frozen, so a caller holding a report line could append a
    pair id to it -- fabricating a changed pair in the record that is supposed
    to be the immutable answer -- and the instance would be unhashable, so
    reports could not be de-duplicated or put in a set. The ids are stored in
    the canonical form `load_pairs` recorded them in, because they are the
    citation a reader compares against cpicpgx.org.
    """

    subject_id: str
    gene: str
    changed_pair_ids: tuple[str, ...]


def affected_by_guideline_change(
    store: RecordStore,
    changed_pair_ids: set[str],
    pairs: list[GuidelineRef],
) -> list[AffectedRecord]:
    """Find stored records involving gene-drug pairs whose guidance changed.

    Returns one entry per (subject, gene), ordered by gene and then by subject
    -- a stable order across runs, independent of the order of `pairs` and of
    the iteration order of `changed_pair_ids`, so two runs of the same revision
    can be diffed. A subject with several stored records for one gene appears
    once: the store is append-only, so a re-ingested gene has many rows and one
    subject.
    """
    # Checked before anything else, and loudly, because the failure is uniform:
    # with no pairs to match against, every revision in existence comes back
    # "nobody affected". `load_pairs` refuses to return an empty table, so this
    # only catches a caller that assembled the list by hand -- but this report
    # is consulted precisely when someone wants to know what moved, and a
    # confident empty answer is the one shape they cannot check.
    if not pairs:
        raise ValueError(
            "refusing to diff against an empty gene-drug pair table: every "
            "guideline revision would be reported as affecting nobody"
        )
    # There is deliberately no `if not changed_pair_ids: return []` short-circuit
    # here. It would be unreachable as a behaviour -- an empty set matches no
    # pair, so the code below already returns [] -- and it would be actively
    # harmful for anything else falsy: `changed_pair_ids=None` from a caller
    # whose "what changed" lookup came back empty-handed would be answered
    # "nobody affected" instead of raising. A caller bug must not be laundered
    # into a confident negative in this particular report.
    #
    # Both sides are normalized, and that is the difference between a report and
    # an empty one. The stored id is canonicalized on load (stripped, case
    # preserved because it is displayed); the caller's id has been through
    # nothing at all, so it arrives however a release note, a spreadsheet cell
    # or a shell argument spelled it. Compared raw, "cyp2c19-clopidogrel" or
    # "CYP2C19-clopidogrel " matches no row, and the function reports that
    # nobody's results changed meaning -- a false negative in the exact report
    # someone opened to find out otherwise. `normalize_pair_id` is the same
    # comparison form `load_pairs` keys its uniqueness check on, so one logical
    # pair cannot be two things here and one thing there.
    changed = {normalize_pair_id(pair_id) for pair_id in changed_pair_ids}

    genes_to_pairs: dict[str, set[str]] = {}
    for pair in pairs:
        if normalize_pair_id(pair.cpic_pair_id) in changed:
            # A gene symbol that is not already canonical cannot match anything:
            # `subjects_with_gene` compares `gene` with SQL `=`, which is
            # case-sensitive for ASCII, and both sides of this system deal in
            # HGNC uppercase (`load_pairs` uppercases the table; PharmCAT emits
            # uppercase symbols, which is what `GeneCall.gene` carries). So a
            # row spelling it "cyp2c19" or " CYP2C19 " would silently contribute
            # zero subjects. Uppercasing it here instead would be worse than
            # raising: it would paper over a table this module cannot see the
            # rest of, while `query_drug` -- which looks the same gene up with an
            # exact dict lookup -- kept answering cannot_assess for the same
            # subjects. Fail loudly, once, where the malformed row is visible.
            if pair.gene != normalize_gene(pair.gene):
                raise ValueError(
                    f"gene-drug pair {pair.cpic_pair_id!r} carries a "
                    f"non-canonical gene symbol {pair.gene!r}; stored calls use "
                    f"HGNC uppercase, so this pair would report no affected "
                    f"subjects even for subjects genotyped for "
                    f"{normalize_gene(pair.gene)!r}"
                )
            # A set, so a table listing one logical pair twice cannot make one
            # changed pair look like two in the report.
            genes_to_pairs.setdefault(pair.gene, set()).add(pair.cpic_pair_id)

    affected: list[AffectedRecord] = []
    for gene, pair_ids in sorted(genes_to_pairs.items()):
        # `subjects_with_gene` already returns distinct subjects, sorted, across
        # every record rather than only the latest -- which is what makes a
        # re-ingested subject one line here and a subject whose only genotyping
        # predates the revision still a line at all.
        for subject_id in store.subjects_with_gene(gene):
            affected.append(
                AffectedRecord(
                    subject_id=subject_id,
                    gene=gene,
                    changed_pair_ids=tuple(sorted(pair_ids)),
                )
            )
    return affected
