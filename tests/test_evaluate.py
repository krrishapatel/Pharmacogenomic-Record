"""Drug query tests.

The invariant under test: absence of guidance and absence of data are
different answers and must never collapse into one. Most of this file exists
to prove that a `cannot_assess` answer cannot be mistaken for an all-clear.

Every path here is anchored on `Path(__file__)`; scratch databases live in
pytest's `tmp_path`. Nothing is read relative to the CWD.
"""

import json
import re
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest

from pharmacogenomic_record import PHARMCAT_VERSION
from pharmacogenomic_record.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall
from pharmacogenomic_record.evaluate import (
    CANNOT_ASSESS,
    GUIDANCE_FOUND,
    NO_GUIDANCE_FOR_PAIR,
    QueryResult,
    _calls_by_gene,
    overall_outcome,
    query_drug,
)
from pharmacogenomic_record.guidelines import (
    GuidelineRef,
    GuidelineTableError,
    find_pairs_for_drug,
    load_pairs,
    normalize_drug,
    normalize_gene,
)
from pharmacogenomic_record.store import CorruptRecordError, RecordStore

PAIRS_PATH = Path(__file__).resolve().parents[1] / "data/gene_drug_pairs.json"
PAIRS = load_pairs(PAIRS_PATH)

GUIDELINE_FIELDS = {"gene", "drug", "cpic_pair_id", "url"}

# Phrases that would turn "we do not know" into an all-clear. Taken from the
# brief and applied to every cannot_assess explanation, not just one.
REASSURING = (
    "no interaction",
    "no issue",
    "safe",
    "normal",
    "no guidance",
    "clear",
    "fine",
)

# Anything that would make the output a recommendation rather than a citation.
DIRECTIVES = ("mg", "you should", "take ", "avoid ", "dose of", "recommend")

# A subject id is caller-supplied text, so it can itself contain every phrase
# the reassurance check looks for. If any explanation interpolates the id, a
# consumer scanning for an all-clear gets one from data the subject named.
HOSTILE_SUBJECT_ID = "patient A - no interaction, safe, normal"

SUBJECT_IDS = [
    pytest.param("s1", id="plain"),
    pytest.param(HOSTILE_SUBJECT_ID, id="hostile"),
]


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


class StubCallSource:
    """A call source that is not a RecordStore, for combinations it forbids.

    `RecordStore` has a CHECK constraint that forbids a phenotype on an uncalled
    row, so through the store an uncalled gene always has phenotype None. That
    makes a guard written as "uncalled AND no phenotype" indistinguishable from
    the correct "uncalled" -- until some other call source (a caller result not
    yet appended, an importer, a future backend) supplies the combination. Only
    `coverage` may decide assessability, so the rule is pinned here.
    """

    def __init__(self, calls):
        self._calls = list(calls)

    def latest(self, subject_id):
        return list(self._calls)


def write_table(tmp_path, payload) -> Path:
    path = tmp_path / "pairs.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


VALID_ENTRY = {
    "gene": "CYP2C19",
    "drug": "clopidogrel",
    "cpic_pair_id": "CYP2C19-clopidogrel",
    "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
}


# --------------------------------------------------------------------------
# The brief's nine tests, verbatim.
# --------------------------------------------------------------------------


def test_guidance_found_for_called_gene(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "clopidogrel", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "guidance_found"
    assert results[0].phenotype == "Intermediate Metabolizer"
    assert results[0].guideline.cpic_pair_id == "CYP2C19-clopidogrel"
    assert results[0].guideline.url.startswith("https://")


def test_no_guidance_for_drug_with_no_cpic_pair(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "amoxicillin", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "no_guidance_for_pair"
    assert results[0].guideline is None


def test_cannot_assess_when_gene_not_covered(store):
    """The invariant. A gene the array never informed carries NO information."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "cannot_assess"
    assert results[0].phenotype is None


def test_cannot_assess_when_gene_absent_from_record(store):
    """A gene missing entirely is unassessable, never 'no interaction'."""
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)

    assert [r.outcome for r in results] == ["cannot_assess"]


def test_cannot_assess_when_indeterminate(store):
    store.append(
        "s1",
        [GeneCall("TPMT", None, None, "indeterminate")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "azathioprine", PAIRS)

    assert results[0].outcome == "cannot_assess"


def test_cannot_assess_is_never_worded_as_reassurance(store):
    """No 'cannot_assess' explanation may read as an all-clear."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)
    text = results[0].explanation.lower()

    for reassuring in (
        "no interaction",
        "no issue",
        "safe",
        "normal",
        "no guidance",
        "clear",
        "fine",
    ):
        assert reassuring not in text, f"reassuring phrase {reassuring!r} in {text!r}"
    assert "cannot" in text or "unable" in text


def test_no_output_contains_a_dose_or_clinical_directive(store):
    """No clinical claims. Guidance is referenced, never restated as advice."""
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "clopidogrel", PAIRS)
    text = results[0].explanation.lower()

    for directive in ("mg", "you should", "take ", "avoid ", "dose of", "recommend"):
        assert directive not in text, f"clinical directive {directive!r} in {text!r}"


def test_unknown_subject_cannot_be_assessed(store):
    results = query_drug(store, "nobody", "clopidogrel", PAIRS)
    assert results[0].outcome == "cannot_assess"


def test_every_outcome_is_one_of_three_values(store):
    store.append(
        "s1",
        [
            GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
            GeneCall("CYP2D6", None, None, "not_covered"),
        ],
        guideline_version="cpic-2026-07",
    )

    outcomes = {
        r.outcome
        for drug in ("clopidogrel", "codeine", "amoxicillin")
        for r in query_drug(store, "s1", drug, PAIRS)
    }
    assert outcomes <= {"guidance_found", "no_guidance_for_pair", "cannot_assess"}


# --------------------------------------------------------------------------
# The outcome vocabulary is a wire contract (the CLI branches on it).
# --------------------------------------------------------------------------


def test_outcome_constants_have_the_documented_values():
    assert GUIDANCE_FOUND == "guidance_found"
    assert NO_GUIDANCE_FOR_PAIR == "no_guidance_for_pair"
    assert CANNOT_ASSESS == "cannot_assess"


# --------------------------------------------------------------------------
# Coverage dominates guidance existence, in every uncalled state.
# --------------------------------------------------------------------------


UNCALLED = [
    pytest.param(NOT_COVERED, id="not_covered"),
    pytest.param(INDETERMINATE, id="indeterminate"),
]


@pytest.mark.parametrize("coverage", UNCALLED)
def test_uncalled_gene_never_yields_guidance_or_a_clean_negative(store, coverage):
    """A drug that HAS CPIC guidance must still come back unassessable."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, coverage)],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)

    # The equality above already excludes both other outcomes; restating it as
    # two != assertions adds nothing.
    assert [r.outcome for r in results] == [CANNOT_ASSESS]
    assert results[0].phenotype is None


@pytest.mark.parametrize("coverage", UNCALLED)
def test_coverage_alone_decides_assessability_even_with_a_phenotype(coverage):
    """`coverage` is the only input to assessability. Never phenotype truthiness.

    A guard written as `coverage != CALLED and phenotype is None` is
    indistinguishable from the correct one through `RecordStore`, whose CHECK
    constraint forbids a phenotype on an uncalled row. Fed from any other call
    source it reports guidance_found for a gene the array never covered -- the
    exact collapse this module forbids. And it cannot be repaired by testing
    phenotype instead: a *called* gene may legitimately have phenotype None
    (F2, F5, VKORC1, CFTR, IFNL3, ABCG2 have no metabolizer phenotype), which
    the null-phenotype test above pins from the other direction.
    """
    source = StubCallSource(
        [GeneCall("CYP2D6", "*1/*1", "Normal Metabolizer", coverage)]
    )

    results = query_drug(source, "s1", "codeine", PAIRS)

    assert [r.outcome for r in results] == [CANNOT_ASSESS]
    # A phenotype on an unassessable gene is not reported: there is no phenotype
    # to stand behind, whatever the row happened to carry.
    assert results[0].phenotype is None
    assert "Normal Metabolizer" not in results[0].explanation


def cannot_assess_situations(tmp_path, subject_id):
    """Every distinct way a query can end in cannot_assess, for one subject id.

    Returned as (label, results) so the reassurance and directive tests can
    cover all of them rather than the single not_covered case in the brief.
    Each situation gets its own store so that every one of them -- including
    the no-record path, which needs a subject the store has never seen -- can
    use the *same* subject id. That is what lets the caller feed a hostile id
    through both the no-record and the record-exists wordings.
    """

    def stored(name, calls):
        store = RecordStore(tmp_path / f"{name}.db")
        if calls is not None:
            store.append(subject_id, calls, guideline_version="cpic-2026-07")
        return store

    covered = stored(
        "covered", [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)]
    )
    uncovered = stored("uncovered", [GeneCall("CYP2D6", None, None, NOT_COVERED)])
    ambiguous = stored("ambiguous", [GeneCall("CYP2D6", None, None, INDETERMINATE)])
    # No append at all: the subject has no stored record of any kind.
    empty = stored("empty", None)
    return [
        ("not_covered", query_drug(uncovered, subject_id, "codeine", PAIRS)),
        ("indeterminate", query_drug(ambiguous, subject_id, "codeine", PAIRS)),
        ("gene absent", query_drug(covered, subject_id, "codeine", PAIRS)),
        ("no record", query_drug(empty, subject_id, "codeine", PAIRS)),
    ]


@pytest.mark.parametrize("subject_id", SUBJECT_IDS)
def test_no_cannot_assess_explanation_reads_as_reassurance(tmp_path, subject_id):
    """The invariant's last line of defense, across all four situations.

    Run twice: once with an ordinary subject id, once with an id that itself
    spells out every reassuring phrase. No explanation may interpolate the
    subject id, because a caller-supplied string inside a safety-critical
    sentence lets the data decide how the answer reads.
    """
    situations = cannot_assess_situations(tmp_path, subject_id)
    assert len(situations) == 4

    for label, results in situations:
        assert [r.outcome for r in results] == [CANNOT_ASSESS], label
        text = results[0].explanation.lower()
        for phrase in REASSURING:
            assert phrase not in text, f"{label}: reassuring {phrase!r} in {text!r}"
        assert "cannot" in text or "unable" in text, f"{label}: {text!r}"
        # It must say, in words, that this is missing data rather than a
        # negative finding.
        assert "unknown" in text or "absence of data" in text, f"{label}: {text!r}"


@pytest.mark.parametrize("subject_id", SUBJECT_IDS)
def test_no_cannot_assess_explanation_quotes_the_subject_id(tmp_path, subject_id):
    """The subject id belongs in the caller's context, not in the explanation.

    Stronger than the reassurance check and independent of it: even an innocuous
    id must stay out, so the leak cannot come back through an id that happens
    not to contain a reassuring word today.
    """
    situations = cannot_assess_situations(tmp_path, subject_id)
    for label, results in situations:
        text = results[0].explanation
        assert subject_id not in text, f"{label}: subject id leaked into {text!r}"
    # Both wordings must still be reachable, or the assertion above is vacuous:
    # "never genotyped" and "genotyped but this gene was not covered" are
    # different states and must not collapse now that the id is gone.
    by_label = dict(situations)
    no_record = by_label["no record"][0].explanation.lower()
    gene_absent = by_label["gene absent"][0].explanation.lower()
    assert no_record != gene_absent
    assert "no record" in no_record
    assert "no call for this gene" in gene_absent


def test_no_explanation_of_any_outcome_contains_a_directive(store):
    """Extends the brief's directive test to every outcome and every drug."""
    store.append(
        "s1",
        [
            GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED),
            GeneCall("CYP2C9", "*1/*3", "Intermediate Metabolizer", CALLED),
            GeneCall("VKORC1", "rs9923231 variant (T)/rs9923231 variant (T)",
                     None, CALLED),
            GeneCall("CYP2D6", None, None, NOT_COVERED),
            GeneCall("TPMT", None, None, INDETERMINATE),
        ],
        guideline_version="cpic-2026-07",
    )

    drugs = sorted({p.drug for p in PAIRS}) + ["amoxicillin", "ibuprofen"]
    seen = set()
    for drug in drugs:
        for result in query_drug(store, "s1", drug, PAIRS):
            seen.add(result.outcome)
            text = result.explanation.lower()
            for directive in DIRECTIVES:
                assert directive not in text, (
                    f"{drug}/{result.outcome}: directive {directive!r} in {text!r}"
                )
    # Only meaningful if all three outcomes were actually exercised.
    assert seen == {GUIDANCE_FOUND, NO_GUIDANCE_FOR_PAIR, CANNOT_ASSESS}


def test_cannot_assess_is_not_interchangeable_with_no_guidance(store):
    """The two negative-looking answers must be distinguishable, always."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )

    unassessable = query_drug(store, "s1", "codeine", PAIRS)[0]
    no_pair = query_drug(store, "s1", "amoxicillin", PAIRS)[0]

    assert unassessable != no_pair
    assert unassessable.outcome != no_pair.outcome
    assert unassessable.explanation != no_pair.explanation
    assert str(unassessable) != str(no_pair)
    # Same outcome field is the only thing a consumer is allowed to branch on,
    # so it must not be possible to make them equal by ignoring the wording.
    assert replace(unassessable, explanation="") != replace(no_pair, explanation="")
    # The unassessable answer names the gene it could not assess; the missing
    # pair answer has no gene to name.
    assert unassessable.gene == "CYP2D6"
    assert no_pair.gene is None


# --------------------------------------------------------------------------
# The citation is the deliverable: every explanation that names a pair must
# carry that pair's URL in its own text.
# --------------------------------------------------------------------------


def test_explanations_cite_the_guideline_url_in_their_own_text(store):
    """Checking only the `guideline.url` field lets the prose drop the link.

    A consumer that renders `explanation` -- which is what the text is for --
    would then show a pair id with nothing to look up, for both the outcome
    that found guidance and the outcome that could not assess it.
    """
    store.append(
        "s1",
        [
            GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED),
            GeneCall("CYP2D6", None, None, NOT_COVERED),
        ],
        guideline_version="cpic-2026-07",
    )

    found = query_drug(store, "s1", "clopidogrel", PAIRS)[0]
    unassessable = query_drug(store, "s1", "codeine", PAIRS)[0]

    assert found.outcome == GUIDANCE_FOUND
    assert found.guideline.url in found.explanation

    # cannot_assess cites the pair too: the reader is told what they *would*
    # consult once the gene is genotyped, not left with a dead end.
    assert unassessable.outcome == CANNOT_ASSESS
    assert unassessable.guideline.url in unassessable.explanation


# --------------------------------------------------------------------------
# A called gene with no metabolizer phenotype is still assessable.
# --------------------------------------------------------------------------


def test_called_gene_with_null_phenotype_is_still_assessable(store):
    """F2/F5/VKORC1/CFTR/IFNL3/ABCG2 have no metabolizer phenotype at all."""
    store.append(
        "s1",
        [
            GeneCall("VKORC1", "rs9923231 variant (T)/rs9923231 variant (T)",
                     None, CALLED),
            GeneCall("CYP2C9", "*1/*3", "Intermediate Metabolizer", CALLED),
        ],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "warfarin", PAIRS)
    by_gene = {r.gene: r for r in results}

    assert by_gene["VKORC1"].outcome == GUIDANCE_FOUND
    assert by_gene["VKORC1"].phenotype is None
    assert by_gene["VKORC1"].guideline.cpic_pair_id == "VKORC1-warfarin"
    assert by_gene["CYP2C9"].outcome == GUIDANCE_FOUND


# --------------------------------------------------------------------------
# Multi-gene drugs: the answer can never be more confident than its
# least-assessable component.
# --------------------------------------------------------------------------


def test_multi_gene_drug_keeps_the_uncovered_gene_visible(store):
    store.append(
        "s1",
        [
            GeneCall("CYP2C9", "*1/*3", "Intermediate Metabolizer", CALLED),
            GeneCall("VKORC1", None, None, NOT_COVERED),
        ],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "warfarin", PAIRS)
    by_gene = {r.gene: r for r in results}

    assert len(results) == 2
    assert by_gene["CYP2C9"].outcome == GUIDANCE_FOUND
    assert by_gene["VKORC1"].outcome == CANNOT_ASSESS
    # Not a clean guidance_found answer: the gap has to survive summarizing.
    assert {r.outcome for r in results} != {GUIDANCE_FOUND}
    assert overall_outcome(results) == CANNOT_ASSESS


def test_overall_outcome_is_the_least_reassuring_component():
    def result(outcome):
        return QueryResult(outcome, "G", None, None, "x")

    assert overall_outcome([result(GUIDANCE_FOUND)]) == GUIDANCE_FOUND
    assert overall_outcome([result(NO_GUIDANCE_FOR_PAIR)]) == NO_GUIDANCE_FOR_PAIR
    assert (
        overall_outcome([result(GUIDANCE_FOUND), result(CANNOT_ASSESS)])
        == CANNOT_ASSESS
    )
    assert (
        overall_outcome([result(NO_GUIDANCE_FOR_PAIR), result(CANNOT_ASSESS)])
        == CANNOT_ASSESS
    )
    # Guidance outranks a missing pair: a real pair must not be hidden behind
    # another gene that CPIC happens not to publish for.
    assert (
        overall_outcome([result(NO_GUIDANCE_FOR_PAIR), result(GUIDANCE_FOUND)])
        == GUIDANCE_FOUND
    )


def test_overall_outcome_refuses_an_empty_or_unknown_answer():
    with pytest.raises(ValueError, match="empty"):
        overall_outcome([])
    with pytest.raises(ValueError, match="unrecognized outcome"):
        overall_outcome([QueryResult("probably_fine", "G", None, None, "x")])


# --------------------------------------------------------------------------
# query_drug never answers with silence.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "drug",
    [
        # One drug in the table and one absent from it -- the two branches of
        # the never-empty guarantee. The trimmed cases were not duplicates of
        # these: with only CYP2C19 stored, "codeine" returns [cannot_assess] and
        # "warfarin" returns two results, both different from "clopidogrel"'s
        # guidance_found. They are covered by the tests named for them --
        # test_cannot_assess_when_gene_absent_from_record and
        # test_multi_gene_drug_keeps_the_uncovered_gene_visible -- which is why
        # dropping them here lost nothing.
        pytest.param("clopidogrel", id="in-table"),
        pytest.param("amoxicillin", id="not-in-table"),
    ],
)
def test_query_drug_never_returns_an_empty_list(store, drug):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    results = query_drug(store, "s1", drug, PAIRS)
    assert results
    assert all(r.outcome in {GUIDANCE_FOUND, NO_GUIDANCE_FOR_PAIR, CANNOT_ASSESS}
               for r in results)


def test_unknown_drug_is_an_explicit_negative_not_an_empty_success(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "not-a-real-drug", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == NO_GUIDANCE_FOR_PAIR
    assert results[0].guideline is None
    assert results[0].gene is None
    # It must say the table is a subset, so "not in our table" is never read
    # as "CPIC publishes nothing".
    text = results[0].explanation.lower()
    assert "not-a-real-drug" in text
    assert "subset" in text or "not the complete" in text


# --------------------------------------------------------------------------
# Distinguishing "no record at all" from "gene absent from a record".
# --------------------------------------------------------------------------


def test_no_record_and_absent_gene_are_both_unassessable_but_distinguishable(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    absent_gene = query_drug(store, "s1", "codeine", PAIRS)[0]
    no_record = query_drug(store, "nobody", "codeine", PAIRS)[0]

    assert absent_gene.outcome == no_record.outcome == CANNOT_ASSESS
    assert absent_gene.explanation != no_record.explanation
    assert "no record" in no_record.explanation.lower()


# --------------------------------------------------------------------------
# Corruption is louder than a coverage gap.
# --------------------------------------------------------------------------


def test_corrupt_record_raises_instead_of_degrading_to_cannot_assess(store):
    """A record with no gene calls is a database integrity failure, not a gap."""
    conn = sqlite3.connect(store.db_path)
    try:
        with conn:
            conn.execute(
                "INSERT INTO records "
                "(subject_id, pharmcat_version, guideline_version, ingested_at) "
                "VALUES (?, ?, ?, ?)",
                ("s1", PHARMCAT_VERSION, "cpic-2026-07",
                 "2026-01-01T00:00:00+00:00"),
            )
    finally:
        conn.close()

    with pytest.raises(CorruptRecordError):
        query_drug(store, "s1", "clopidogrel", PAIRS)


# --------------------------------------------------------------------------
# Drug name matching.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "drug",
    ["clopidogrel", "Clopidogrel", "CLOPIDOGREL", "  clopidogrel  ", "\tClopidogreL\n"],
)
def test_drug_matching_ignores_case_and_surrounding_whitespace(store, drug):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", drug, PAIRS)

    assert [r.outcome for r in results] == [GUIDANCE_FOUND]
    assert results[0].guideline.cpic_pair_id == "CYP2C19-clopidogrel"


@pytest.mark.parametrize(
    "drug",
    [
        # Only these two actually require NFKC. Nothing else in the codebase
        # folds a fullwidth letter, so removing the normalize() call breaks
        # exactly these params -- which is what makes them worth keeping.
        pytest.param("ＷＡＲＦＡＲＩＮ", id="fullwidth-upper"),
        pytest.param("ｗａｒｆａｒｉｎ", id="fullwidth-lower"),
        # Kept, but honestly labelled: `"\xa0".strip()` is already "" in
        # Python, so this param passes with or without NFKC. It pins that
        # stripping covers non-ASCII whitespace; it demonstrates nothing about
        # NFKC, and its old id ("non-breaking-space") implied otherwise.
        pytest.param(" warfarin ", id="nbsp-padding-stripped-without-nfkc"),
        # What a real paste actually looks like: fullwidth letters *inside*
        # NBSP padding. Needs NFKC for the letters and strip() for the
        # padding, so it fails if either half is removed.
        pytest.param(
            " ｗａｒｆａｒｉｎ ",
            id="fullwidth-inside-nbsp-padding",
        ),
    ],
)
def test_compatibility_forms_of_a_drug_name_still_match(store, drug):
    """A name pasted from a PDF or typed on an IME must not answer 'no guidance'.

    Falling through to no_guidance_for_pair is the safe direction, but it is
    still a wrong answer the user has no reason to question.
    """
    store.append(
        "s1",
        [
            GeneCall("CYP2C9", "*1/*3", "Intermediate Metabolizer", CALLED),
            GeneCall("VKORC1", "rs9923231 variant (T)/rs9923231 variant (T)",
                     None, CALLED),
        ],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", drug, PAIRS)

    assert {r.gene for r in results} == {"CYP2C9", "VKORC1"}
    assert {r.outcome for r in results} == {GUIDANCE_FOUND}


def test_normalization_does_not_collide_two_distinct_shipped_drugs():
    """NFKC must fold compatibility forms without merging different drugs."""
    drugs = {p.drug for p in PAIRS}
    normalized = {normalize_drug(d) for d in drugs}
    assert len(normalized) == len(drugs)
    # Each shipped name is already in normal form, so normalizing is a no-op.
    for drug in drugs:
        assert normalize_drug(drug) == drug


def test_a_homoglyph_of_a_drug_name_fails_closed():
    """Cyrillic "а" is not Latin "a" and NFKC does not fold it -- correctly.

    A spoofed name must not match a real drug. Refusing to recognize the name
    is the right failure: the alternative is matching the wrong drug.
    """
    spoofed = "wаrfarin"
    assert spoofed != "warfarin"
    assert normalize_drug(spoofed) != normalize_drug("warfarin")
    assert find_pairs_for_drug(spoofed, PAIRS) == []


@pytest.mark.parametrize("drug", ["", "   ", "\t\n", " ", "　"])
def test_blank_drug_name_is_refused_not_answered(store, drug):
    """A blank query answered 'no guidance' is a fabricated negative."""
    with pytest.raises(ValueError, match="blank"):
        query_drug(store, "s1", drug, PAIRS)


def test_non_string_drug_name_is_refused(store):
    with pytest.raises(TypeError):
        query_drug(store, "s1", None, PAIRS)


def test_empty_pair_table_is_refused_not_answered(store):
    """An empty table would report every drug on earth as uncovered by CPIC."""
    with pytest.raises(ValueError, match="empty gene-drug pair table"):
        query_drug(store, "s1", "clopidogrel", [])


def test_find_pairs_for_drug_returns_every_gene_for_a_multi_gene_drug():
    genes = {p.gene for p in find_pairs_for_drug("WARFARIN", PAIRS)}
    assert genes == {"CYP2C9", "VKORC1"}
    assert find_pairs_for_drug("amoxicillin", PAIRS) == []


def test_gene_matching_errs_toward_cannot_assess(store):
    """Gene names are matched exactly; a mismatch fails safe, never confident."""
    store.append(
        "s1",
        [GeneCall("cyp2c19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "clopidogrel", PAIRS)

    assert [r.outcome for r in results] == [CANNOT_ASSESS]


def test_duplicate_gene_with_conflicting_coverage_is_refused():
    """Two rows for one gene must not silently resolve to the confident one.

    Unreachable through RecordStore -- (record_id, gene) is a primary key -- so
    the helper is exercised directly. The check exists because a dict built by
    last-write-wins would let a `called` row mask a `not_covered` one.
    """
    calls = [
        GeneCall("CYP2D6", None, None, NOT_COVERED),
        GeneCall("CYP2D6", "*1/*1", "Normal Metabolizer", CALLED),
    ]
    with pytest.raises(CorruptRecordError, match="CYP2D6"):
        _calls_by_gene(calls)

    assert _calls_by_gene(calls[:1]) == {"CYP2D6": calls[0]}


# --------------------------------------------------------------------------
# The shipped pair table: identifiers and links only, never prose.
# --------------------------------------------------------------------------


def test_shipped_table_stores_identifiers_and_links_only():
    entries = json.loads(PAIRS_PATH.read_text(encoding="utf-8"))
    assert isinstance(entries, list)
    assert entries

    for entry in entries:
        assert set(entry) == GUIDELINE_FIELDS, entry
        assert entry["url"].startswith("https://cpicpgx.org/"), entry
        for key, value in entry.items():
            assert isinstance(value, str) and value.strip(), entry
            # A prose field would be long; identifiers and URLs are not.
            assert len(value) < 200, key


def test_shipped_table_covers_the_pairs_the_cli_and_drift_tests_use():
    ids = {p.cpic_pair_id for p in PAIRS}
    assert {
        "CYP2C19-clopidogrel",
        "CYP2D6-codeine",
        "DPYD-fluorouracil",
        "TPMT-azathioprine",
        "SLCO1B1-simvastatin",
        "CYP2C9-warfarin",
        "VKORC1-warfarin",
    } <= ids
    assert len(ids) == len(PAIRS)


def test_load_pairs_round_trips_a_minimal_table(tmp_path):
    path = write_table(tmp_path, [VALID_ENTRY])
    assert load_pairs(path) == [GuidelineRef(**VALID_ENTRY)]


@pytest.mark.parametrize(
    "payload, match",
    [
        pytest.param({}, "list", id="object-not-list"),
        pytest.param([], "no gene-drug pairs", id="empty-list"),
        pytest.param(["CYP2C19-clopidogrel"], "object", id="string-entry"),
        pytest.param(
            [{k: v for k, v in VALID_ENTRY.items() if k != "url"}],
            "url",
            id="missing-field",
        ),
        pytest.param(
            [{**VALID_ENTRY, "recommendation": "reduce the dose"}],
            "unexpected",
            id="extra-prose-field",
        ),
        pytest.param([{**VALID_ENTRY, "gene": "  "}], "blank", id="blank-field"),
        pytest.param([{**VALID_ENTRY, "drug": 7}], "string", id="non-string-field"),
        pytest.param(
            [{**VALID_ENTRY, "url": "http://cpicpgx.org/x"}],
            "https",
            id="insecure-url",
        ),
        pytest.param([VALID_ENTRY, VALID_ENTRY], "duplicate", id="duplicate-pair"),
        # The licensing guard. It exists so nobody pastes guideline prose into a
        # field, which is the thing the module docstring says this table never
        # carries -- and which we have no verified right to redistribute.
        pytest.param(
            [{**VALID_ENTRY, "cpic_pair_id": "CYP2C19-clopidogrel " + "x" * 200}],
            "not guideline prose",
            id="over-long-field",
        ),
        # A citation that does not match its own row: the wrong guideline, cited
        # confidently, is a wrong answer nothing downstream can detect.
        pytest.param(
            [{**VALID_ENTRY, "cpic_pair_id": "SLCO1B1-simvastatin"}],
            "does not name its own gene",
            id="pair-id-for-another-gene",
        ),
        pytest.param(
            [{**VALID_ENTRY, "url": "https://example.com/x"}],
            "cpicpgx.org",
            id="url-off-host",
        ),
        # A host that merely *contains* the real one must not pass.
        pytest.param(
            [{**VALID_ENTRY, "url": "https://cpicpgx.org.evil.example/x"}],
            "cpicpgx.org",
            id="url-lookalike-host",
        ),
    ],
)
def test_load_pairs_rejects_an_untrustworthy_table(tmp_path, payload, match):
    """A silently truncated table turns 'guidance exists' into 'no guidance'."""
    path = write_table(tmp_path, payload)
    with pytest.raises(GuidelineTableError, match=match):
        load_pairs(path)


def test_a_typod_citation_is_refused_rather_than_answered(tmp_path):
    """The reviewer's exact row: right gene and drug, someone else's guideline.

    Loaded, this answers guidance_found for CYP2C19/clopidogrel while citing the
    statin guideline at an arbitrary host. The whole output of this tool is a
    citation, so a wrong citation is a wrong answer that looks entirely
    confident. Nothing loads at all instead.
    """
    path = write_table(
        tmp_path,
        [
            {
                "gene": "CYP2C19",
                "drug": "clopidogrel",
                "cpic_pair_id": "SLCO1B1-simvastatin",
                "url": "https://example.com/x",
            }
        ],
    )
    with pytest.raises(GuidelineTableError, match="CYP2C19"):
        load_pairs(path)


def test_the_shipped_table_still_loads_under_every_citation_check():
    """All seven rows satisfy the gene/pair-id and host checks as shipped."""
    assert len(PAIRS) == 7
    for pair in PAIRS:
        assert pair.gene.casefold() in pair.cpic_pair_id.casefold(), pair
        # Not merely a substring: the gene must sit on a token boundary in its
        # own pair id, which is the check load_pairs actually applies.
        assert re.search(
            rf"(?<![0-9A-Za-z]){re.escape(pair.gene)}(?![0-9A-Za-z])",
            pair.cpic_pair_id,
            re.IGNORECASE,
        ), pair
        assert pair.url.startswith("https://cpicpgx.org/"), pair


# --------------------------------------------------------------------------
# The gene/pair-id check matches a token, not a bare substring.
# --------------------------------------------------------------------------


def test_a_gene_embedded_in_a_longer_symbol_is_not_a_citation_match(tmp_path):
    """"F2" is a substring of "CYP4F2-warfarin", and that is not a match.

    F2 and CYP4F2 are both real warfarin-associated genes, and F2 is in this
    project's own null-phenotype list, so this row is a plausible typo rather
    than a contrived one. Under a bare substring test it LOADED, and a called F2
    then came back guidance_found citing the CYP4F2 guideline -- the exact silent
    wrong-gene citation this check exists to stop, and one nothing downstream can
    catch, because the citation IS the answer.
    """
    path = write_table(
        tmp_path,
        [
            {
                "gene": "F2",
                "drug": "warfarin",
                "cpic_pair_id": "CYP4F2-warfarin",
                "url": "https://cpicpgx.org/guidelines/"
                "guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
            }
        ],
    )
    with pytest.raises(GuidelineTableError, match="does not name its own gene"):
        load_pairs(path)


@pytest.mark.parametrize(
    "gene, cpic_pair_id",
    [
        pytest.param("CYP2C19", "CYP2C19-clopidogrel", id="gene-first"),
        pytest.param("CYP2C9", "CYP2C9-warfarin", id="shorter-cyp2c-sibling"),
        pytest.param("VKORC1", "VKORC1-warfarin", id="trailing-digit"),
        pytest.param(
            "DPYD",
            "guideline-for-fluoropyrimidines-and-dpyd",
            id="gene-last-and-lowercase",
        ),
        pytest.param("HLA-B", "HLA-B-abacavir", id="hyphen-inside-the-symbol"),
        pytest.param("F2", "F2-warfarin", id="two-character-symbol"),
        pytest.param("SLCO1B1", "SLCO1B1-simvastatin", id="digits-inside"),
        pytest.param("TPMT", "TPMT-azathioprine", id="letters-only"),
        pytest.param("CYP2D6", "CYP2D6-codeine", id="cyp2d6"),
    ],
)
def test_every_legitimate_pair_id_shape_still_loads(tmp_path, gene, cpic_pair_id):
    """The accept side of the boundary check, pinned shape by shape.

    A rejected row is not a rejected row: load_pairs is all-or-nothing, so one
    of these failing means the whole table stops loading and EVERY drug answers
    no_guidance_for_pair -- a fabricated negative across the board. This list is
    what stops a future "tightening" of the check from doing that silently.
    """
    path = write_table(
        tmp_path,
        [{**VALID_ENTRY, "gene": gene, "cpic_pair_id": cpic_pair_id}],
    )
    loaded = load_pairs(path)
    assert [p.cpic_pair_id for p in loaded] == [cpic_pair_id]
    assert loaded[0].gene == gene


# --------------------------------------------------------------------------
# A gene symbol is stored in the form every comparison will use.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "gene",
    [
        pytest.param(" CYP2C19 ", id="padded"),
        pytest.param("cyp2c19", id="lowercase"),
        pytest.param("\tCyp2c19\n", id="padded-and-mixed-case"),
    ],
)
def test_a_padded_or_case_variant_gene_is_canonicalized_not_stored_as_written(
    tmp_path, store, gene
):
    """Storing the raw value made the row validate and then never match.

    The old check compared `gene.strip()` but stored the unstripped string, so a
    " CYP2C19 " row loaded happily and a genuinely called CYP2C19 came back
    cannot_assess -- `query_drug` looks the gene up by exact dict key. The safe
    direction, but still a confident wrong answer about a subject who WAS
    genotyped, and the strip inside the check is what hid it.
    """
    path = write_table(tmp_path, [{**VALID_ENTRY, "gene": gene}])
    pairs = load_pairs(path)
    assert [p.gene for p in pairs] == ["CYP2C19"]

    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    results = query_drug(store, "s1", "clopidogrel", pairs)
    assert [r.outcome for r in results] == [GUIDANCE_FOUND]


@pytest.mark.parametrize(
    "second_gene",
    [
        pytest.param("cyp2c19", id="case-variant"),
        pytest.param(" CYP2C19 ", id="padded"),
    ],
)
def test_a_case_variant_gene_cannot_evade_the_duplicate_pair_check(
    tmp_path, second_gene
):
    """Uniqueness is keyed on the gene, so it has to be keyed on ONE form.

    With the raw value stored, "CYP2C19" and "cyp2c19" were different keys and a
    second copy of an existing pair slipped through -- reported twice for one
    query. The pair ids here are deliberately distinct so only the gene key can
    catch it.
    """
    path = write_table(
        tmp_path,
        [
            VALID_ENTRY,
            {
                **VALID_ENTRY,
                "gene": second_gene,
                "cpic_pair_id": "CYP2C19-clopidogrel-2019",
            },
        ],
    )
    with pytest.raises(GuidelineTableError, match="duplicate gene-drug pair"):
        load_pairs(path)


def test_normalize_gene_is_the_canonical_form_used_on_load():
    assert normalize_gene(" cyp2c19 ") == "CYP2C19"
    assert normalize_gene("HLA-B") == "HLA-B"
    # Every shipped gene is already canonical, so loading is a no-op for them.
    for pair in PAIRS:
        assert normalize_gene(pair.gene) == pair.gene


# --------------------------------------------------------------------------
# A url this table cannot parse is still named as one bad row.
# --------------------------------------------------------------------------


def test_an_unparseable_url_is_reported_as_a_table_error_naming_the_entry(tmp_path):
    """`urlsplit` raises a bare ValueError, which escapes the error contract.

    "https://[abc/x" raises ValueError("Invalid IPv6 URL") with no entry index
    and no mention of the pair table, so the one failure a maintainer cannot fix
    from the message is the one that says least. Every other rejection in
    load_pairs names its row.
    """
    path = write_table(
        tmp_path, [VALID_ENTRY, {**VALID_ENTRY, "url": "https://[abc/x"}]
    )
    with pytest.raises(GuidelineTableError) as excinfo:
        load_pairs(path)
    message = str(excinfo.value)
    assert "entry 1" in message, message
    assert "https://[abc/x" in message, message


def test_only_the_gene_is_pinned_to_the_pair_id_not_the_drug(tmp_path):
    """CPIC names some guidelines after a drug class rather than a member.

    "DPYD-fluoropyrimidines" is a legitimate pair id for a row whose drug is
    capecitabine, so requiring the drug to appear in the pair id would reject a
    correct row -- and a rejected row loads as nothing at all, which is the
    false negative this module exists to avoid.
    """
    path = write_table(
        tmp_path,
        [
            {
                "gene": "DPYD",
                "drug": "capecitabine",
                "cpic_pair_id": "DPYD-fluoropyrimidines",
                "url": "https://cpicpgx.org/guidelines/"
                "guideline-for-fluoropyrimidines-and-dpyd/",
            }
        ],
    )
    assert [p.drug for p in load_pairs(path)] == ["capecitabine"]


def test_load_pairs_rejects_unreadable_and_unparseable_files(tmp_path):
    with pytest.raises(GuidelineTableError, match="could not read"):
        load_pairs(tmp_path / "missing.json")

    broken = tmp_path / "broken.json"
    broken.write_text("[{", encoding="utf-8")
    with pytest.raises(GuidelineTableError, match="could not read"):
        load_pairs(broken)


def test_guideline_ref_is_frozen():
    ref = GuidelineRef(**VALID_ENTRY)
    with pytest.raises(Exception):
        ref.url = "https://example.com/"


def test_query_result_is_frozen(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    result = query_drug(store, "s1", "clopidogrel", PAIRS)[0]
    with pytest.raises(Exception):
        result.outcome = CANNOT_ASSESS
