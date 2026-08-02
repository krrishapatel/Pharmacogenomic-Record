import sqlite3
from datetime import datetime, timedelta, timezone

import pytest

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall
from pgxrecord.store import RecordStore

CALLS = [
    GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
    GeneCall("CYP2D6", None, None, "not_covered"),
]


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


def test_append_and_read_back(store):
    store.append("subject-1", CALLS, guideline_version="cpic-2026-07")
    calls = store.latest("subject-1")

    by_gene = {c.gene: c for c in calls}
    assert by_gene["CYP2C19"].phenotype == "Intermediate Metabolizer"
    assert by_gene["CYP2D6"].coverage == "not_covered"


def test_versions_are_stamped(store):
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    pharmcat_version, guideline_version = store.record_versions(record_id)

    assert pharmcat_version == PHARMCAT_VERSION
    assert guideline_version == "cpic-2026-07"


def test_reingest_appends_and_preserves_history(store):
    first = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    updated = [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called")]
    second = store.append("s1", updated, guideline_version="cpic-2026-08")

    assert store.history("s1") == [first, second]
    assert store.latest("s1")[0].phenotype == "Normal Metabolizer"
    # The original call is still retrievable -- history is not destroyed.
    assert store.record_versions(first)[1] == "cpic-2026-07"


def test_updates_are_rejected_at_the_database_level(store):
    """Immutability is enforced by the schema, not by convention."""
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("UPDATE gene_calls SET phenotype = 'tampered'")
        conn.commit()
    conn.close()


def test_deletes_are_rejected_at_the_database_level(store):
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("DELETE FROM gene_calls")
        conn.commit()
    conn.close()


def test_subjects_with_gene_finds_affected_records(store):
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    store.append("s2", [GeneCall("DPYD", "Ref/Ref", "Normal", "called")],
                 guideline_version="cpic-2026-07")

    assert store.subjects_with_gene("CYP2C19") == ["s1"]
    assert store.subjects_with_gene("DPYD") == ["s2"]


def test_unknown_subject_returns_empty_history(store):
    assert store.history("nobody") == []
    assert store.latest("nobody") == []


# --- Immutability, in more detail than the brief's two trigger tests ---------


def test_rejected_update_leaves_the_stored_call_untouched(store):
    """A rejected UPDATE must not have partially rewritten anything."""
    store.append("s1", CALLS, guideline_version="cpic-2026-07")

    conn = sqlite3.connect(store.db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("UPDATE gene_calls SET phenotype = 'tampered'")
    conn.close()

    by_gene = {c.gene: c for c in store.latest("s1")}
    assert by_gene["CYP2C19"].phenotype == "Intermediate Metabolizer"
    assert by_gene["CYP2D6"].phenotype is None


def test_records_table_rejects_update_and_delete(store):
    """The version stamp is as immutable as the calls it explains."""
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute("UPDATE records SET guideline_version = 'cpic-1999-01'")
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute("UPDATE records SET pharmcat_version = '0.0.0'")
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute("DELETE FROM records")
    finally:
        conn.close()

    assert store.record_versions(store.history("s1")[0]) == (
        PHARMCAT_VERSION,
        "cpic-2026-07",
    )


def test_insert_or_replace_cannot_overwrite_a_stored_call(store):
    """REPLACE is an UPDATE wearing an INSERT's clothes, and must also fail.

    sqlite implements INSERT OR REPLACE as a delete-then-insert, but it skips
    the BEFORE DELETE trigger unless `PRAGMA recursive_triggers` is on -- and
    that pragma is per-connection and off by default. A caller who opens the
    file directly therefore gets the unsafe setting for free, so append-only
    cannot rest on the delete trigger alone.
    """
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        for statement in (
            "INSERT OR REPLACE INTO gene_calls "
            "(record_id, gene, diplotype, phenotype, coverage) "
            "VALUES (?, ?, ?, ?, ?)",
            "REPLACE INTO gene_calls "
            "(record_id, gene, diplotype, phenotype, coverage) "
            "VALUES (?, ?, ?, ?, ?)",
        ):
            with pytest.raises(sqlite3.IntegrityError):
                conn.execute(
                    statement,
                    (record_id, "CYP2C19", "*17/*17", "Ultrarapid Metabolizer",
                     CALLED),
                )
                conn.commit()
    finally:
        conn.close()

    by_gene = {c.gene: c for c in store.record_calls(record_id)}
    assert by_gene["CYP2C19"].phenotype == "Intermediate Metabolizer"
    assert by_gene["CYP2C19"].diplotype == "*1/*2"


def test_insert_or_replace_cannot_overwrite_a_version_stamp(store):
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT OR REPLACE INTO records (record_id, subject_id, "
                "pharmcat_version, guideline_version, ingested_at) "
                "VALUES (?, ?, ?, ?, ?)",
                (record_id, "s1", "0.0.0", "tampered",
                 "2020-01-01T00:00:00+00:00"),
            )
            conn.commit()
    finally:
        conn.close()

    assert store.record_versions(record_id) == (PHARMCAT_VERSION, "cpic-2026-07")


def test_appending_a_second_call_for_a_gene_in_the_same_record_is_rejected(store):
    """A record is a snapshot: one call per gene, and it is closed once written."""
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO gene_calls (record_id, gene, diplotype, phenotype, "
                "coverage) VALUES (?, ?, ?, ?, ?)",
                (record_id, "CYP2C19", "*2/*2", "Poor Metabolizer", CALLED),
            )
            conn.commit()
    finally:
        conn.close()

    assert len(store.record_calls(record_id)) == len(CALLS)


def test_reingest_keeps_both_entries_retrievable_and_the_earlier_unchanged(store):
    """Two entries, both readable, the older one bit-for-bit as written."""
    first = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    updated = [
        GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called"),
        GeneCall("CYP2D6", "*1/*4", "Intermediate Metabolizer", "called"),
    ]
    second = store.append("s1", updated, guideline_version="cpic-2026-08")

    assert first != second
    assert store.record_calls(first) == sorted(CALLS, key=lambda c: c.gene)
    assert store.record_calls(second) == sorted(updated, key=lambda c: c.gene)
    assert store.record_versions(first) == (PHARMCAT_VERSION, "cpic-2026-07")
    assert store.record_versions(second) == (PHARMCAT_VERSION, "cpic-2026-08")
    # The formerly-uncovered gene is now called, and the record that said
    # "we did not know" is still there to explain the change.
    assert {c.gene: c.coverage for c in store.record_calls(first)}["CYP2D6"] == (
        NOT_COVERED
    )


def test_history_is_oldest_to_newest_across_many_appends(store):
    ids = [
        store.append("s1", CALLS, guideline_version=f"cpic-2026-{month:02d}")
        for month in range(1, 6)
    ]
    assert store.history("s1") == ids
    assert ids == sorted(ids)


# --- Atomicity: a failed write leaves no partial entry -----------------------


def test_failed_write_leaves_no_partial_entry(store):
    """The gene_calls insert fails after the records insert; both roll back."""
    bad = [
        GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
        GeneCall("CYP2D6", None, None, "probably_fine"),
    ]
    with pytest.raises(sqlite3.IntegrityError):
        store.append("s1", bad, guideline_version="cpic-2026-07")

    assert store.history("s1") == []
    assert store.latest("s1") == []
    assert store.subjects_with_gene("CYP2C19") == []


def test_failed_write_does_not_disturb_earlier_entries(store):
    good = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    bad = [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "guesswork")]
    with pytest.raises(sqlite3.IntegrityError):
        store.append("s1", bad, guideline_version="cpic-2026-08")

    assert store.history("s1") == [good]
    assert store.latest("s1")[0].phenotype == "Intermediate Metabolizer"


def test_duplicate_gene_in_one_record_is_rejected_atomically(store):
    """One record, one call per gene -- and a rejected record is no record."""
    dupes = [
        GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
        GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called"),
    ]
    with pytest.raises(sqlite3.IntegrityError):
        store.append("s1", dupes, guideline_version="cpic-2026-07")

    assert store.history("s1") == []


def test_empty_call_list_is_rejected(store):
    """An entry with no calls is unreadable: latest() would look like 'no data'."""
    with pytest.raises(ValueError):
        store.append("s1", [], guideline_version="cpic-2026-07")

    assert store.history("s1") == []


# --- Version stamping is structural, not incidental -------------------------


def test_every_stored_entry_has_non_null_versions(store):
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    store.append("s2", CALLS, guideline_version="cpic-2026-08")

    conn = sqlite3.connect(store.db_path)
    try:
        missing = conn.execute(
            "SELECT COUNT(*) FROM records WHERE pharmcat_version IS NULL "
            "OR guideline_version IS NULL OR trim(pharmcat_version) = '' "
            "OR trim(guideline_version) = ''"
        ).fetchone()[0]
        total = conn.execute("SELECT COUNT(*) FROM records").fetchone()[0]
    finally:
        conn.close()

    assert total == 2
    assert missing == 0


def test_unstamped_record_cannot_be_inserted_at_all(store):
    """The schema, not the application, is what makes the stamp mandatory."""
    conn = sqlite3.connect(store.db_path)
    try:
        for pharmcat_version, guideline_version in [
            (None, "cpic-2026-07"),
            (PHARMCAT_VERSION, None),
            ("", "cpic-2026-07"),
            (PHARMCAT_VERSION, ""),
            (PHARMCAT_VERSION, "   "),
        ]:
            with pytest.raises(sqlite3.IntegrityError):
                conn.execute(
                    "INSERT INTO records (subject_id, pharmcat_version, "
                    "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
                    ("s1", pharmcat_version, guideline_version, "2026-07-31T00:00:00+00:00"),
                )
    finally:
        conn.close()


def test_blank_guideline_version_is_rejected_by_append(store):
    with pytest.raises(sqlite3.IntegrityError):
        store.append("s1", CALLS, guideline_version="  ")
    assert store.history("s1") == []


# --- Coverage states survive the round trip ---------------------------------


def test_all_three_coverage_states_round_trip_distinctly(store):
    calls = [
        GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED),
        GeneCall("CYP2D6", None, None, NOT_COVERED),
        GeneCall("NAT2", None, None, INDETERMINATE),
    ]
    store.append("s1", calls, guideline_version="cpic-2026-07")

    by_gene = {c.gene: c for c in store.latest("s1")}
    assert by_gene["CYP2C19"].coverage == CALLED
    assert by_gene["CYP2D6"].coverage == NOT_COVERED
    assert by_gene["NAT2"].coverage == INDETERMINATE
    assert len({c.coverage for c in by_gene.values()}) == 3

    # not_covered and indeterminate must not collapse into each other, and
    # neither may come back carrying a call.
    assert by_gene["CYP2D6"].coverage != by_gene["NAT2"].coverage
    for gene in ("CYP2D6", "NAT2"):
        assert by_gene[gene].diplotype is None
        assert by_gene[gene].phenotype is None

    assert by_gene["CYP2C19"] == calls[0]


def test_a_gene_absent_from_a_record_is_absent_not_called(store):
    """`latest()` reports only genes the record mentions -- and no default.

    Downstream code must read a missing gene as "we have no record for this
    gene", never as a negative finding. There is deliberately no filled-in
    placeholder row that could be mistaken for a result.
    """
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    genes = {c.gene for c in store.latest("s1")}
    assert genes == {"CYP2C19"}
    assert "CYP2D6" not in genes


def test_called_gene_with_no_phenotype_round_trips(store):
    """F2/F5/VKORC1 have no metabolizer phenotype; that is not indeterminate."""
    store.append(
        "s1",
        [GeneCall("VKORC1", "rs9923231 variant (T)/rs9923231 variant (T)", None, CALLED)],
        guideline_version="cpic-2026-07",
    )
    call = store.latest("s1")[0]
    assert call.coverage == CALLED
    assert call.phenotype is None
    assert call.diplotype == "rs9923231 variant (T)/rs9923231 variant (T)"


def test_unknown_coverage_state_is_rejected(store):
    with pytest.raises(sqlite3.IntegrityError):
        store.append(
            "s1",
            [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "probable")],
            guideline_version="cpic-2026-07",
        )


def test_coverage_cannot_be_smuggled_past_the_check_by_type_coercion(store):
    """sqlite's dynamic typing must not open a path to a bogus coverage value.

    A TEXT column happily stores an integer or a blob, so the CHECK -- not the
    declared type -- is what keeps the coverage vocabulary closed.
    """
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        for bogus in (None, 1, 1.0, b"called", "CALLED", " called", ""):
            with pytest.raises(sqlite3.IntegrityError):
                conn.execute(
                    "INSERT INTO gene_calls (record_id, gene, diplotype, "
                    "phenotype, coverage) VALUES (?, ?, ?, ?, ?)",
                    (record_id, "GBA", None, None, bogus),
                )
    finally:
        conn.close()


def test_uncalled_gene_cannot_carry_a_call_even_via_raw_sql(store):
    """not_covered/indeterminate rows are structurally forbidden a diplotype."""
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    try:
        for diplotype, phenotype, coverage in [
            ("*1/*1", None, NOT_COVERED),
            (None, "Normal Metabolizer", NOT_COVERED),
            ("*1/*1", "Normal Metabolizer", INDETERMINATE),
            (None, "Normal Metabolizer", INDETERMINATE),
        ]:
            with pytest.raises(sqlite3.IntegrityError):
                conn.execute(
                    "INSERT INTO gene_calls (record_id, gene, diplotype, "
                    "phenotype, coverage) VALUES (?, ?, ?, ?, ?)",
                    (record_id, "GBA", diplotype, phenotype, coverage),
                )
    finally:
        conn.close()


def test_called_gene_without_a_diplotype_is_rejected(store):
    """A "called" row with nothing to report is a contradiction, not a call."""
    for diplotype in (None, "", "   "):
        with pytest.raises(sqlite3.IntegrityError):
            store.append(
                "s1",
                [GeneCall("CYP2C19", diplotype, "Normal Metabolizer", CALLED)],
                guideline_version="cpic-2026-07",
            )
    assert store.history("s1") == []


# --- Timestamps -------------------------------------------------------------


def test_ingested_at_is_injectable_and_read_back_exactly(store):
    when = datetime(2026, 7, 31, 12, 34, 56, tzinfo=timezone.utc)
    record_id = store.append(
        "s1", CALLS, guideline_version="cpic-2026-07", ingested_at=when
    )
    assert store.record_ingested_at(record_id) == when


def test_ingested_at_defaults_to_now_in_utc(store):
    before = datetime.now(timezone.utc)
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    after = datetime.now(timezone.utc)

    stamped = store.record_ingested_at(record_id)
    assert stamped.tzinfo is not None
    assert before - timedelta(seconds=1) <= stamped <= after + timedelta(seconds=1)


def test_naive_ingested_at_is_rejected(store):
    with pytest.raises(ValueError):
        store.append(
            "s1",
            CALLS,
            guideline_version="cpic-2026-07",
            ingested_at=datetime(2026, 7, 31, 12, 0, 0),
        )
    assert store.history("s1") == []


def test_ingested_at_is_normalized_to_utc(store):
    offset = timezone(timedelta(hours=-7))
    when = datetime(2026, 7, 31, 5, 34, 56, tzinfo=offset)
    record_id = store.append(
        "s1", CALLS, guideline_version="cpic-2026-07", ingested_at=when
    )
    stamped = store.record_ingested_at(record_id)
    assert stamped == when
    assert stamped.utcoffset() == timedelta(0)


# --- Lookups ----------------------------------------------------------------


def test_unknown_record_id_raises(store):
    with pytest.raises(KeyError):
        store.record_versions(999)
    with pytest.raises(KeyError):
        store.record_calls(999)
    with pytest.raises(KeyError):
        store.record_ingested_at(999)


def test_subjects_with_gene_spans_all_records_not_just_the_latest(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )
    store.append(
        "s1",
        [GeneCall("DPYD", "Ref/Ref", "Normal Metabolizer", CALLED)],
        guideline_version="cpic-2026-08",
    )
    # A guideline change for CYP2C19 still concerns s1, even though the most
    # recent record does not mention the gene.
    assert store.subjects_with_gene("CYP2C19") == ["s1"]
    assert store.subjects_with_gene("nonesuch") == []


def test_subjects_with_gene_deduplicates_and_sorts(store):
    for subject in ("s2", "s1", "s2"):
        store.append(subject, CALLS, guideline_version="cpic-2026-07")
    assert store.subjects_with_gene("CYP2D6") == ["s1", "s2"]


def test_a_store_reopened_on_the_same_path_sees_the_same_history(store, tmp_path):
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    reopened = RecordStore(tmp_path / "records.db")

    assert reopened.history("s1") == [record_id]
    assert reopened.record_calls(record_id) == sorted(CALLS, key=lambda c: c.gene)
    assert reopened.record_versions(record_id)[1] == "cpic-2026-07"


def test_store_does_not_leak_connections_across_many_appends(tmp_path):
    """Every helper must close its handle; sqlite3's own context manager does not."""
    store = RecordStore(tmp_path / "records.db")
    for month in range(1, 60):
        store.append("s1", CALLS, guideline_version=f"g-{month}")
        store.latest("s1")
        store.history("s1")
    assert len(store.history("s1")) == 59
