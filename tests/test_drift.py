"""Guideline-drift tests.

The invariant under test: when guidance moves, every subject the change could
possibly concern must appear in the report -- including subjects whose gene was
never covered, because "we never looked at your CYP2D6" is precisely the reason
new CYP2D6 guidance might warrant testing. A report that quietly omits someone
is worse than no report, because its emptiness reads as reassurance.

Two silent-zero failure modes get their own tests: a caller-supplied pair id
that differs only in case or padding, and a coverage state other than `called`.
Both would return an empty list -- indistinguishable from "nobody affected".

Every path here is anchored on `Path(__file__)`; scratch databases live in
pytest's `tmp_path`. Nothing is read relative to the CWD.
"""

import dataclasses
from pathlib import Path

import pytest

from pharmacogenomic_record.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall
from pharmacogenomic_record.drift import AffectedRecord, affected_by_guideline_change
from pharmacogenomic_record.guidelines import GuidelineRef, load_pairs
from pharmacogenomic_record.store import RecordStore

PAIRS_PATH = Path(__file__).resolve().parents[1] / "data/gene_drug_pairs.json"
PAIRS = load_pairs(PAIRS_PATH)

VERSION = "cpic-2026-07"


@pytest.fixture
def store(tmp_path):
    store = RecordStore(tmp_path / "records.db")
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version=VERSION,
    )
    store.append(
        "s2",
        [GeneCall("DPYD", "Ref/Ref", "Normal Metabolizer", CALLED)],
        guideline_version=VERSION,
    )
    store.append(
        "s3",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version=VERSION,
    )
    return store


# --- the brief's six ---------------------------------------------------------


def test_finds_subjects_affected_by_a_changed_pair(store):
    affected = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)

    assert len(affected) == 1
    assert affected[0].subject_id == "s1"
    assert affected[0].gene == "CYP2C19"
    assert affected[0].changed_pair_ids == ("CYP2C19-clopidogrel",)


def test_unaffected_subjects_are_not_reported(store):
    affected = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)
    assert {a.subject_id for a in affected} == {"s1"}


def test_multiple_changed_pairs(store):
    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel", "DPYD-fluorouracil"}, PAIRS
    )
    assert {a.subject_id for a in affected} == {"s1", "s2"}


def test_not_covered_genes_are_still_reported(store):
    """A not_covered gene still matters: new guidance may warrant re-testing."""
    affected = affected_by_guideline_change(store, {"CYP2D6-codeine"}, PAIRS)

    assert len(affected) == 1
    assert affected[0].subject_id == "s3"
    assert affected[0].gene == "CYP2D6"


def test_unknown_pair_id_affects_nobody(store):
    assert affected_by_guideline_change(store, {"NOPE-nothing"}, PAIRS) == []


def test_no_changes_affects_nobody(store):
    assert affected_by_guideline_change(store, set(), PAIRS) == []


# --- the silent zero: a caller-supplied id that is not byte-identical --------


@pytest.mark.parametrize(
    "supplied",
    [
        pytest.param("cyp2c19-clopidogrel", id="lowercase"),
        pytest.param("CYP2C19-CLOPIDOGREL", id="uppercase"),
        pytest.param("CYP2C19-Clopidogrel", id="titlecase"),
        pytest.param(" CYP2C19-clopidogrel ", id="padded"),
        pytest.param("\tCYP2C19-clopidogrel\n", id="tab-and-newline"),
        pytest.param(" cyp2c19-CLOPIDOGREL\t", id="padded-and-miscased"),
    ],
)
def test_caller_pair_id_case_and_padding_still_finds_the_subject(store, supplied):
    """A miscased or padded id must not report "nobody affected".

    The stored id is canonicalized (stripped) on load while the caller's id has
    been through nothing, so a raw `in` test against the stored value misses.
    Missing here does not raise -- it returns [], which reads exactly like a
    revision that touched none of these records.
    """
    affected = affected_by_guideline_change(store, {supplied}, PAIRS)

    assert [(a.subject_id, a.gene) for a in affected] == [("s1", "CYP2C19")]


def test_reported_ids_are_the_stored_form_not_the_callers(store):
    """The id we report is the citation a user checks against cpicpgx.org."""
    affected = affected_by_guideline_change(store, {"cyp2c19-clopidogrel"}, PAIRS)

    assert affected[0].changed_pair_ids == ("CYP2C19-clopidogrel",)


def test_case_variant_ids_collapse_rather_than_duplicating(store):
    """Two spellings of one id are one changed pair, not two."""
    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel", "cyp2c19-clopidogrel "}, PAIRS
    )

    assert len(affected) == 1
    assert affected[0].changed_pair_ids == ("CYP2C19-clopidogrel",)


# --- coverage is never a filter ---------------------------------------------


@pytest.mark.parametrize(
    "coverage",
    [
        pytest.param(NOT_COVERED, id="not_covered"),
        pytest.param(INDETERMINATE, id="indeterminate"),
    ],
)
def test_uncalled_coverage_states_are_all_reported(tmp_path, coverage):
    """Absence of data is not absence of an interaction, in either flavour."""
    store = RecordStore(tmp_path / f"{coverage}.db")
    store.append(
        "s1", [GeneCall("CYP2C19", None, None, coverage)], guideline_version=VERSION
    )

    affected = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)

    assert [(a.subject_id, a.gene) for a in affected] == [("s1", "CYP2C19")]


def test_a_called_and_an_uncalled_subject_are_reported_together(tmp_path):
    """The uncalled subject must not be dropped from a report that has hits."""
    store = RecordStore(tmp_path / "mixed.db")
    store.append(
        "called",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version=VERSION,
    )
    store.append(
        "uncalled",
        [GeneCall("CYP2C19", None, None, NOT_COVERED)],
        guideline_version=VERSION,
    )

    affected = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)

    assert [a.subject_id for a in affected] == ["called", "uncalled"]


# --- determinism -------------------------------------------------------------


@pytest.fixture
def crowd(tmp_path):
    """Several subjects across several genes, appended out of sorted order.

    The genes span the shipped table widely enough that the order the report
    groups them in -- the pair table's order -- is not their sorted order. The
    table lists CYP2C19, DPYD, TPMT, SLCO1B1, CYP2C9; sorted they run CYP2C19,
    CYP2C9, DPYD, SLCO1B1, TPMT. A report that skipped the gene sort would
    therefore come out visibly reordered rather than passing by coincidence.
    """
    store = RecordStore(tmp_path / "crowd.db")
    for subject_id, gene in [
        ("s2", "DPYD"),
        ("s1", "CYP2C19"),
        ("s3", "CYP2C19"),
        ("s2", "CYP2C19"),
        ("s1", "DPYD"),
        ("s3", "DPYD"),
        ("s2", "TPMT"),
        ("s3", "SLCO1B1"),
        ("s1", "CYP2C9"),
        ("s2", "SLCO1B1"),
        ("s3", "CYP2C9"),
        ("s1", "TPMT"),
        ("s3", "TPMT"),
        ("s1", "SLCO1B1"),
        ("s2", "CYP2C9"),
    ]:
        store.append(
            subject_id,
            [GeneCall(gene, "*1/*1", "Normal Metabolizer", CALLED)],
            guideline_version=VERSION,
        )
    return store


def test_order_is_gene_then_subject_and_repeats_identically(crowd):
    """Genes sorted, then subjects sorted -- not the pair table's own order.

    The changed pairs are supplied for five genes whose table order (CYP2C19,
    DPYD, TPMT, SLCO1B1, CYP2C9) is not their sorted order, so dropping the gene
    sort and grouping in table-encounter order fails this rather than
    coincidentally satisfying it.
    """
    changed = {
        "DPYD-fluorouracil",
        "CYP2C19-clopidogrel",
        "TPMT-azathioprine",
        "SLCO1B1-simvastatin",
        "CYP2C9-warfarin",
    }

    affected = affected_by_guideline_change(crowd, changed, PAIRS)

    assert [(a.gene, a.subject_id) for a in affected] == [
        ("CYP2C19", "s1"),
        ("CYP2C19", "s2"),
        ("CYP2C19", "s3"),
        ("CYP2C9", "s1"),
        ("CYP2C9", "s2"),
        ("CYP2C9", "s3"),
        ("DPYD", "s1"),
        ("DPYD", "s2"),
        ("DPYD", "s3"),
        ("SLCO1B1", "s1"),
        ("SLCO1B1", "s2"),
        ("SLCO1B1", "s3"),
        ("TPMT", "s1"),
        ("TPMT", "s2"),
        ("TPMT", "s3"),
    ]
    # Same store, same inputs, same list -- twice, and independent of the order
    # the caller's set happens to iterate in.
    assert affected_by_guideline_change(crowd, changed, PAIRS) == affected
    assert (
        affected_by_guideline_change(crowd, set(reversed(sorted(changed))), PAIRS)
        == affected
    )


def test_order_is_independent_of_the_pair_table_order(crowd):
    changed = {"CYP2C19-clopidogrel", "DPYD-fluorouracil"}

    forward = affected_by_guideline_change(crowd, changed, PAIRS)
    backward = affected_by_guideline_change(crowd, changed, list(reversed(PAIRS)))

    assert forward == backward


# --- one drug, two genes; one gene, two pairs -------------------------------


def test_subject_with_two_changed_genes_is_reported_once_per_gene(tmp_path):
    """Warfarin is two rows: CYP2C9-warfarin and VKORC1-warfarin."""
    store = RecordStore(tmp_path / "warfarin.db")
    store.append(
        "s1",
        [
            GeneCall("CYP2C9", "*1/*2", "Intermediate Metabolizer", CALLED),
            GeneCall("VKORC1", "rs9923231 variant", None, CALLED),
        ],
        guideline_version=VERSION,
    )

    affected = affected_by_guideline_change(
        store, {"CYP2C9-warfarin", "VKORC1-warfarin"}, PAIRS
    )

    assert affected == [
        AffectedRecord("s1", "CYP2C9", ("CYP2C9-warfarin",)),
        AffectedRecord("s1", "VKORC1", ("VKORC1-warfarin",)),
    ]


VORICONAZOLE = GuidelineRef(
    gene="CYP2C19",
    drug="voriconazole",
    cpic_pair_id="CYP2C19-voriconazole",
    url="https://cpicpgx.org/guidelines/guideline-for-voriconazole-and-cyp2c19/",
)


@pytest.mark.parametrize(
    "pairs",
    [
        pytest.param(PAIRS + [VORICONAZOLE], id="table-order-matches-sorted"),
        pytest.param([VORICONAZOLE] + PAIRS, id="table-order-reversed"),
    ],
)
def test_one_gene_in_two_changed_pairs_groups_both_ids(store, pairs):
    """Both of a gene's changed pairs land on its single AffectedRecord.

    Parametrized on where the second pair sits in the table, because the ids
    must come back sorted rather than in whatever order the table happens to
    list them: a report whose lines reshuffle between runs cannot be diffed.
    """
    affected = affected_by_guideline_change(
        store, {"CYP2C19-voriconazole", "CYP2C19-clopidogrel"}, pairs
    )

    assert affected == [
        AffectedRecord(
            "s1", "CYP2C19", ("CYP2C19-clopidogrel", "CYP2C19-voriconazole")
        )
    ]


def test_grouped_pair_ids_are_sorted_not_in_table_order(store):
    """Enough ids that an unsorted grouping cannot coincidentally come out sorted.

    Two ids collected into a set land in sorted order for a good fraction of
    hash seeds, so a two-id assertion lets an unsorted implementation pass on
    those runs. These five are supplied in reverse and asserted forward.
    """
    extra = [
        GuidelineRef(
            gene="CYP2C19",
            drug=drug,
            cpic_pair_id=f"CYP2C19-{drug}",
            url=f"https://cpicpgx.org/guidelines/guideline-for-{drug}-and-cyp2c19/",
        )
        for drug in ["voriconazole", "sertraline", "escitalopram", "lansoprazole"]
    ]
    changed = {"CYP2C19-clopidogrel"} | {p.cpic_pair_id for p in extra}

    affected = affected_by_guideline_change(store, changed, PAIRS + extra)

    ids = affected[0].changed_pair_ids
    assert ids == tuple(sorted(ids))
    assert ids == (
        "CYP2C19-clopidogrel",
        "CYP2C19-escitalopram",
        "CYP2C19-lansoprazole",
        "CYP2C19-sertraline",
        "CYP2C19-voriconazole",
    )


def test_only_the_changed_pair_of_a_multi_pair_gene_is_reported(store):
    affected = affected_by_guideline_change(
        store, {"CYP2C19-voriconazole"}, PAIRS + [VORICONAZOLE]
    )

    assert affected == [AffectedRecord("s1", "CYP2C19", ("CYP2C19-voriconazole",))]


# --- append-only store, one report line ------------------------------------


def test_repeated_records_for_one_gene_report_the_subject_once(tmp_path):
    """The store is append-only, so one gene can appear in many records."""
    store = RecordStore(tmp_path / "reingested.db")
    for diplotype in ["*1/*2", "*1/*2", "*2/*2"]:
        store.append(
            "s1",
            [GeneCall("CYP2C19", diplotype, "Intermediate Metabolizer", CALLED)],
            guideline_version=VERSION,
        )

    affected = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)

    assert affected == [AffectedRecord("s1", "CYP2C19", ("CYP2C19-clopidogrel",))]


def test_an_empty_store_affects_nobody(tmp_path):
    """Empty because nobody is stored -- and it stops being empty when they are.

    The bare `[]` assertion is a boundary case no mutation can fail, since an
    implementation that always returned `[]` would satisfy it. Appending one
    subject to the same store and re-asserting pins the emptiness to the store's
    contents, and the guards still apply to an empty store: it is not a
    shortcut past them.
    """
    store = RecordStore(tmp_path / "empty.db")

    assert affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS) == []
    with pytest.raises(ValueError, match="empty gene-drug pair table"):
        affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, [])
    with pytest.raises(TypeError, match="not a single str"):
        affected_by_guideline_change(store, "CYP2C19-clopidogrel", PAIRS)

    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version=VERSION,
    )

    assert affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS) == [
        AffectedRecord("s1", "CYP2C19", ("CYP2C19-clopidogrel",))
    ]


# --- the report line itself -------------------------------------------------


def test_affected_record_is_frozen(store):
    record = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)[0]

    with pytest.raises(dataclasses.FrozenInstanceError):
        record.subject_id = "someone else"
    with pytest.raises(dataclasses.FrozenInstanceError):
        record.changed_pair_ids = ("CYP2D6-codeine",)


def test_affected_record_is_hashable_and_its_ids_cannot_be_mutated(store):
    record = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)[0]

    # A list field would make the frozen guarantee shallow (and the record
    # unhashable): a caller holding the report could append a pair id to it.
    assert isinstance(record.changed_pair_ids, tuple)
    assert {record, record} == {record}
    with pytest.raises(AttributeError):
        record.changed_pair_ids.append("CYP2D6-codeine")


def test_report_carries_ids_only_and_no_prose(store):
    """Output is a reference: identifiers, never guidance text or a dose."""
    record = affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, PAIRS)[0]

    assert [f.name for f in dataclasses.fields(record)] == [
        "subject_id",
        "gene",
        "changed_pair_ids",
    ]


# --- an empty table would report nobody affected, uniformly ----------------


def test_empty_pair_table_is_refused(store):
    """Every drug's guidance would look unchanged for everyone. Silently."""
    with pytest.raises(ValueError, match="empty gene-drug pair table"):
        affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, [])


def test_empty_pair_table_is_refused_even_when_nothing_changed(store):
    """No falsy short-circuit may run ahead of the empty-table guard.

    Not a duplicate of the test above: this is the case an
    `if not changed_pair_ids: return []` at the top of the function would break,
    turning a refusal into "nobody affected". The empty change set itself is not
    refused -- with a real table it returns [] -- only the empty table is.
    """
    with pytest.raises(ValueError, match="empty gene-drug pair table"):
        affected_by_guideline_change(store, set(), [])

    assert affected_by_guideline_change(store, set(), PAIRS) == []


def test_none_change_set_is_not_answered_as_nobody_affected(store):
    """A caller whose "what changed" lookup returned None gets an error, not [].

    `None` is not "nothing changed"; it is "we do not know what changed", and
    those must not share an answer in this report.
    """
    with pytest.raises(TypeError):
        affected_by_guideline_change(store, None, PAIRS)


# --- one id passed bare instead of wrapped ----------------------------------


def test_a_bare_string_change_set_is_refused_not_reported_as_nobody(store):
    """`"CYP2C9-warfarin"` instead of `{"CYP2C9-warfarin"}` must not return [].

    Iterating a string yields its characters, so the normalized change set
    becomes a bag of letters that matches no pair and the report comes back
    empty -- for a store that genuinely holds an affected subject. This is the
    likeliest caller mistake and the one failure this module must never have,
    so it raises where `None` already did.
    """
    store.append(
        "s4",
        [GeneCall("CYP2C9", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version=VERSION,
    )

    with pytest.raises(TypeError, match="not a single str"):
        affected_by_guideline_change(store, "CYP2C9-warfarin", PAIRS)

    # The same id, wrapped, is not an empty report -- so the [] a bare string
    # used to produce was a fabricated negative, not the truth.
    assert affected_by_guideline_change(store, {"CYP2C9-warfarin"}, PAIRS) == [
        AffectedRecord("s4", "CYP2C9", ("CYP2C9-warfarin",))
    ]


def test_a_bare_bytes_change_set_is_refused(store):
    """Bytes iterate into ints, which `normalize_pair_id` would blow up on.

    Refused with the same TypeError as `str` rather than an AttributeError from
    somewhere deeper, so the message names what the caller actually did wrong.
    """
    with pytest.raises(TypeError, match="not a single bytes"):
        affected_by_guideline_change(store, b"CYP2C19-clopidogrel", PAIRS)


@pytest.mark.parametrize(
    "wrap",
    [
        pytest.param(set, id="set"),
        pytest.param(frozenset, id="frozenset"),
        pytest.param(list, id="list"),
        pytest.param(tuple, id="tuple"),
        pytest.param(lambda ids: (i for i in ids), id="generator"),
        pytest.param(lambda ids: dict.fromkeys(ids).keys(), id="dict-keys"),
        pytest.param(lambda ids: iter(list(ids)), id="iterator"),
    ],
)
def test_every_container_of_ids_still_works(store, wrap):
    """Rejecting `str` must not narrow what a legitimate caller may pass.

    Any iterable of ids worked before the guard and must still work after it,
    including the one-shot iterables a "what changed" query returns.
    """
    affected = affected_by_guideline_change(store, wrap(["CYP2C19-clopidogrel"]), PAIRS)

    assert affected == [AffectedRecord("s1", "CYP2C19", ("CYP2C19-clopidogrel",))]


# --- a gene symbol that could never match a stored call ----------------------


@pytest.mark.parametrize(
    "gene",
    [
        pytest.param("cyp2c19", id="lowercase"),
        pytest.param(" CYP2C19 ", id="padded"),
    ],
)
def test_non_canonical_gene_symbol_is_refused_not_silently_unmatched(store, gene):
    """`subjects_with_gene` compares with SQL `=`, so this would find nobody.

    `load_pairs` canonicalizes `gene`, so such a row can only arrive from a
    hand-built table -- and it would contribute zero affected subjects while
    looking like a pair that was checked.
    """
    hand_built = [
        GuidelineRef(
            gene=gene,
            drug="clopidogrel",
            cpic_pair_id="CYP2C19-clopidogrel",
            url="https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
        )
    ]

    with pytest.raises(ValueError, match="non-canonical gene symbol"):
        affected_by_guideline_change(store, {"CYP2C19-clopidogrel"}, hand_built)


def test_non_canonical_gene_in_an_unchanged_pair_is_ignored(store):
    """Only pairs this report actually consults are held to the guard."""
    hand_built = PAIRS + [
        GuidelineRef(
            gene="cyp2b6",
            drug="efavirenz",
            cpic_pair_id="CYP2B6-efavirenz",
            url="https://cpicpgx.org/guidelines/cpic-guideline-for-efavirenz-based-regimens/",
        )
    ]

    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel"}, hand_built
    )

    assert [a.subject_id for a in affected] == ["s1"]


def test_the_shipped_pair_table_has_only_canonical_gene_symbols(store):
    """The guard above must be unreachable for the table we actually ship."""
    every_id = {p.cpic_pair_id for p in PAIRS}

    affected = affected_by_guideline_change(store, every_id, PAIRS)

    assert [(a.gene, a.subject_id) for a in affected] == [
        ("CYP2C19", "s1"),
        ("CYP2D6", "s3"),
        ("DPYD", "s2"),
    ]
