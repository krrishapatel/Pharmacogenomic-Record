"""Drug query tests.

The invariant under test: absence of guidance and absence of data are
different answers and must never collapse into one. Most of this file exists
to prove that a `cannot_assess` answer cannot be mistaken for an all-clear.

Every path here is anchored on `Path(__file__)`; scratch databases live in
pytest's `tmp_path`. Nothing is read relative to the CWD.
"""

import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall
from pgxrecord.evaluate import (
    CANNOT_ASSESS,
    GUIDANCE_FOUND,
    NO_GUIDANCE_FOR_PAIR,
    QueryResult,
    _calls_by_gene,
    overall_outcome,
    query_drug,
)
from pgxrecord.guidelines import (
    GuidelineRef,
    GuidelineTableError,
    find_pairs_for_drug,
    load_pairs,
)
from pgxrecord.store import CorruptRecordError, RecordStore

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


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


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

    assert [r.outcome for r in results] == [CANNOT_ASSESS]
    assert results[0].outcome != GUIDANCE_FOUND
    assert results[0].outcome != NO_GUIDANCE_FOR_PAIR
    assert results[0].phenotype is None


def cannot_assess_situations(store):
    """Every distinct way a query can end in cannot_assess.

    Returned as (label, results) so the reassurance and directive tests can
    cover all of them rather than the single not_covered case in the brief.
    """
    store.append(
        "covered",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    store.append(
        "uncovered",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )
    store.append(
        "ambiguous",
        [GeneCall("CYP2D6", None, None, INDETERMINATE)],
        guideline_version="cpic-2026-07",
    )
    return [
        ("not_covered", query_drug(store, "uncovered", "codeine", PAIRS)),
        ("indeterminate", query_drug(store, "ambiguous", "codeine", PAIRS)),
        ("gene absent", query_drug(store, "covered", "codeine", PAIRS)),
        ("no record", query_drug(store, "nobody", "codeine", PAIRS)),
    ]


def test_no_cannot_assess_explanation_reads_as_reassurance(store):
    """The invariant's last line of defense, across all four situations."""
    situations = cannot_assess_situations(store)
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
    "drug", ["clopidogrel", "codeine", "warfarin", "amoxicillin", "zzz-not-a-drug"]
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


@pytest.mark.parametrize("drug", ["", "   ", "\t\n"])
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
    ],
)
def test_load_pairs_rejects_an_untrustworthy_table(tmp_path, payload, match):
    """A silently truncated table turns 'guidance exists' into 'no guidance'."""
    path = write_table(tmp_path, payload)
    with pytest.raises(GuidelineTableError, match=match):
        load_pairs(path)


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
