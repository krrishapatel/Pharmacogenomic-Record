import os
import sqlite3
from datetime import datetime, timedelta, timezone
import pytest

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall
from pgxrecord.store import CorruptRecordError, RecordStore

# Every scratch database lives in pytest's tmp_path; this file reads no
# fixtures from disk, so there is deliberately no fixtures path here. The
# other test modules anchor theirs on Path(__file__), never on the CWD.

CALLS = [
    GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
    GeneCall("CYP2D6", None, None, "not_covered"),
]

# Every character sqlite's two-argument trim() is handed in the schema. The
# one-argument trim() strips SPACES ONLY, so all but the last of these sailed
# through a `trim(x) <> ''` CHECK.
WHITESPACE = [
    pytest.param("\t", id="tab"),
    pytest.param("\n", id="newline"),
    pytest.param("\r", id="carriage-return"),
    pytest.param("\v", id="vertical-tab"),
    pytest.param("\f", id="form-feed"),
    pytest.param(" ", id="space"),
]


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


def raw(store):
    """A connection that has never heard of RecordStore, for tamper tests."""
    return sqlite3.connect(store.db_path)


def insert_bare_record(store, subject_id, ingested_at="2026-01-01T00:00:00+00:00"):
    """Write a `records` row directly, with no gene_calls. Returns its id.

    Unreachable through `append()` -- which is the point: these are the shapes
    a corrupt file can have, and the store must refuse to read them rather than
    answer with something plausible.
    """
    conn = raw(store)
    try:
        cursor = conn.execute(
            "INSERT INTO records (subject_id, pharmcat_version, "
            "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
            (subject_id, PHARMCAT_VERSION, "cpic-2026-07", ingested_at),
        )
        record_id = int(cursor.lastrowid)
        conn.commit()
    finally:
        conn.close()
    return record_id


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


def test_insert_or_replace_on_an_explicit_rowid_cannot_overwrite_a_call(store):
    """The rowid back door: REPLACE targeting a rowid, not the (record_id, gene) key.

    `gene_calls` used to be an ordinary rowid table. `INSERT OR REPLACE` naming
    an explicit rowid resolved the conflict on the *rowid*, so:

    * the internal delete skipped BEFORE DELETE (sqlite only fires it under
      `PRAGMA recursive_triggers`, which is per-connection and off by default);
    * `gene_calls_no_replace` checks for an existing (record_id, gene) and the
      replacement row carried a *different* gene, so it never fired either.

    One statement therefore deleted a stored CYP2C19 phenotype and substituted a
    fabricated gene, past all six triggers. `WITHOUT ROWID` removes the key that
    made this addressable, so the statement can no longer be written at all --
    hence sqlite3.OperationalError ("no column named rowid") rather than an
    IntegrityError. Both are sqlite3.Error; what matters is that it is refused
    and that nothing moved.
    """
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")

    conn = raw(store)
    try:
        # Find the rowid an attacker would aim at. On a WITHOUT ROWID table the
        # column does not exist, which is itself the fix -- so tolerate that and
        # fall back to the value the old schema would have handed out.
        try:
            target_rowid = conn.execute(
                "SELECT rowid FROM gene_calls WHERE record_id = ? AND gene = ?",
                (record_id, "CYP2C19"),
            ).fetchone()[0]
        except sqlite3.OperationalError:
            target_rowid = 1

        with pytest.raises(sqlite3.Error) as excinfo:
            conn.execute(
                "INSERT OR REPLACE INTO gene_calls "
                "(rowid, record_id, gene, diplotype, phenotype, coverage) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (target_rowid, record_id, "ZZZ", "*9/*9", "NM", CALLED),
            )
            conn.commit()
        assert isinstance(excinfo.value, sqlite3.Error)
    finally:
        conn.close()

    # The stored history is bit-for-bit what append() wrote: nothing deleted,
    # nothing fabricated.
    by_gene = {c.gene: c for c in store.record_calls(record_id)}
    assert set(by_gene) == {"CYP2C19", "CYP2D6"}
    assert "ZZZ" not in by_gene
    assert by_gene["CYP2C19"] == CALLS[0]
    assert by_gene["CYP2D6"] == CALLS[1]
    assert store.record_calls(record_id) == sorted(CALLS, key=lambda c: c.gene)
    assert store.subjects_with_gene("ZZZ") == []


def test_gene_calls_has_no_rowid_to_target(store):
    """The absence of the rowid is the mechanism, so assert on it directly.

    A future migration that drops WITHOUT ROWID for "performance" reopens the
    hole above, and would do so without breaking any behavioural test that does
    not name the rowid. This one breaks.
    """
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = raw(store)
    try:
        with pytest.raises(sqlite3.OperationalError):
            conn.execute("SELECT rowid FROM gene_calls").fetchall()
        # `records` keeps its rowid on purpose -- record_id *is* the rowid there,
        # so a REPLACE naming it collides on the primary key and the
        # records_no_replace trigger sees it.
        assert conn.execute("SELECT rowid FROM records").fetchall()
    finally:
        conn.close()


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


def test_history_is_insertion_order_even_when_stamps_are_out_of_order(store):
    """A backdated replay must not rewrite what the store considers most recent.

    Every other ordering test appends in increasing stamp order, so
    `ORDER BY record_id` and `ORDER BY ingested_at` agree and neither is pinned.
    Here the second append is stamped a year EARLIER than the first, which makes
    the two orderings disagree: `history()` must stay in insertion order and
    `latest()` must return the record inserted last, not the one with the newest
    stamp. Importing an old archive is exactly this shape, and it must not
    silently demote a newer result nor promote an older one.
    """
    recent = datetime(2026, 7, 31, 12, 0, 0, tzinfo=timezone.utc)
    backdated = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)

    first = store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
        ingested_at=recent,
    )
    second = store.append(
        "s1",
        [GeneCall("CYP2C19", "*2/*2", "Poor Metabolizer", CALLED)],
        guideline_version="cpic-2024-01",
        ingested_at=backdated,
    )

    # The stamps really are inverted relative to insertion order, otherwise the
    # rest of this test would prove nothing.
    assert store.record_ingested_at(first) > store.record_ingested_at(second)
    assert first < second

    assert store.history("s1") == [first, second]
    # ORDER BY ingested_at would return [second, first] here.
    assert store.history("s1") != [second, first]

    # latest() follows insertion order too, not the clock.
    assert store.latest("s1") == store.record_calls(second)
    assert store.latest("s1")[0].phenotype == "Poor Metabolizer"
    assert store.record_calls(first)[0].phenotype == "Normal Metabolizer"


def test_a_backdated_append_does_not_hide_behind_an_existing_record(store):
    """Three appends, middle one backdated: insertion order is still total."""
    stamps = [
        datetime(2026, 6, 1, tzinfo=timezone.utc),
        datetime(2020, 1, 1, tzinfo=timezone.utc),
        datetime(2023, 3, 3, tzinfo=timezone.utc),
    ]
    ids = [
        store.append(
            "s1", CALLS, guideline_version=f"g-{i}", ingested_at=stamp
        )
        for i, stamp in enumerate(stamps)
    ]

    assert store.history("s1") == ids
    chronological = sorted(ids, key=store.record_ingested_at)
    assert chronological != ids, "stamps must disagree with insertion order"
    assert store.record_versions(store.history("s1")[-1])[1] == "g-2"


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


# --- Blankness means all whitespace, not just spaces ------------------------
#
# sqlite's one-argument trim() strips SPACES ONLY, so `trim(char(9)) <> ''` is
# 1 and a `trim(x) <> ''` CHECK accepts a tab, newline, carriage return,
# vertical tab or form feed as a subject id, a gene name or a guideline
# version. A record whose guideline_version is a tab does not know which
# guidance produced it -- exactly what the constraint exists to prevent -- and
# a gene named "\n" is indistinguishable from a gene named "\t" in any report.


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_guideline_version_is_rejected(store, blank):
    with pytest.raises(sqlite3.IntegrityError):
        store.append("s1", CALLS, guideline_version=blank)
    assert store.history("s1") == []


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_guideline_version_is_rejected_via_raw_sql(store, blank):
    """The CHECK, not append(), is what holds. Assert it from outside the module."""
    conn = raw(store)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO records (subject_id, pharmcat_version, "
                "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
                ("s1", PHARMCAT_VERSION, blank, "2026-07-31T00:00:00+00:00"),
            )
    finally:
        conn.close()
    assert store.history("s1") == []


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_pharmcat_version_is_rejected(store, blank):
    conn = raw(store)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO records (subject_id, pharmcat_version, "
                "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
                ("s1", blank, "cpic-2026-07", "2026-07-31T00:00:00+00:00"),
            )
    finally:
        conn.close()


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_subject_id_is_rejected(store, blank):
    """A record nobody can be identified with is not a record."""
    with pytest.raises(sqlite3.IntegrityError):
        store.append(blank, CALLS, guideline_version="cpic-2026-07")
    assert store.history(blank) == []

    conn = raw(store)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO records (subject_id, pharmcat_version, "
                "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
                (blank, PHARMCAT_VERSION, "cpic-2026-07",
                 "2026-07-31T00:00:00+00:00"),
            )
        assert conn.execute("SELECT COUNT(*) FROM records").fetchone()[0] == 0
    finally:
        conn.close()


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_gene_is_rejected(store, blank):
    """A call attributed to no gene cannot be reported or re-evaluated."""
    with pytest.raises(sqlite3.IntegrityError):
        store.append(
            "s1",
            [GeneCall(blank, None, None, NOT_COVERED)],
            guideline_version="cpic-2026-07",
        )
    assert store.history("s1") == []

    record_id = store.append("s2", CALLS, guideline_version="cpic-2026-07")
    conn = raw(store)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO gene_calls (record_id, gene, diplotype, phenotype, "
                "coverage) VALUES (?, ?, ?, ?, ?)",
                (record_id, blank, None, None, NOT_COVERED),
            )
    finally:
        conn.close()
    assert len(store.record_calls(record_id)) == len(CALLS)


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_diplotype_on_a_called_gene_is_rejected(store, blank):
    """"called" plus a tab for a diplotype is a confident result with no content."""
    with pytest.raises(sqlite3.IntegrityError):
        store.append(
            "s1",
            [GeneCall("CYP2C19", blank, "Normal Metabolizer", CALLED)],
            guideline_version="cpic-2026-07",
        )
    assert store.history("s1") == []


@pytest.mark.parametrize("blank", WHITESPACE)
def test_whitespace_only_ingested_at_is_rejected(store, blank):
    conn = raw(store)
    try:
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO records (subject_id, pharmcat_version, "
                "guideline_version, ingested_at) VALUES (?, ?, ?, ?)",
                ("s1", PHARMCAT_VERSION, "cpic-2026-07", blank),
            )
    finally:
        conn.close()


def test_a_gene_name_of_mixed_whitespace_is_rejected(store):
    """Blankness is not about any single character but about having no content."""
    with pytest.raises(sqlite3.IntegrityError):
        store.append(
            "s1",
            [GeneCall(" \t\r\n\v\f ", None, None, INDETERMINATE)],
            guideline_version="cpic-2026-07",
        )
    assert store.history("s1") == []


def test_a_gene_name_containing_whitespace_is_still_accepted(store):
    """Only *entirely* blank is rejected; internal whitespace is not our business.

    Guards against a fix that over-corrects into stripping or rejecting any
    value with whitespace in it. Real PharmCAT diplotype labels contain spaces
    ("rs9923231 variant (T)"), so a whitespace-phobic CHECK would reject valid
    clinical data.
    """
    record_id = store.append(
        "s1",
        [GeneCall("HLA-B", "*15:02 variant", "Positive", CALLED)],
        guideline_version="cpic 2026 07",
    )
    call = store.record_calls(record_id)[0]
    assert call.gene == "HLA-B"
    assert call.diplotype == "*15:02 variant"
    assert store.record_versions(record_id)[1] == "cpic 2026 07"


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
    """UTC normalization must be asserted on the normalized *form*, not on ==.

    `stamped == when` is true in every zone -- aware datetimes compare on the
    instant -- and `utcoffset() == 0` is true whenever Python parses an offset
    back, even one that is not UTC. Neither notices if `astimezone(utc)` is
    dropped from `append` or from `record_ingested_at`. So this asserts on the
    two things that actually change: the string on disk, and the wall clock
    fields that come back.
    """
    offset = timezone(timedelta(hours=-8))
    when = datetime(2026, 1, 1, 0, 0, 0, tzinfo=offset)
    record_id = store.append(
        "s1", CALLS, guideline_version="cpic-2026-07", ingested_at=when
    )

    # What is written to disk is UTC: a file read by a different tool, or by a
    # SQL string comparison, must not have to know about the writer's zone.
    conn = raw(store)
    try:
        stored = conn.execute(
            "SELECT ingested_at FROM records WHERE record_id = ?", (record_id,)
        ).fetchone()[0]
    finally:
        conn.close()
    assert stored.endswith("+00:00"), stored
    assert stored == "2026-01-01T08:00:00+00:00"

    stamped = store.record_ingested_at(record_id)
    assert stamped == when
    assert stamped.utcoffset() == timedelta(0)
    assert stamped.tzinfo is timezone.utc
    # Midnight in UTC-8 is 08:00 UTC. The reader hands back UTC wall-clock
    # fields, not the submitter's, so anything formatting this stamp reports
    # one canonical time rather than the ingesting machine's local one.
    assert (stamped.year, stamped.month, stamped.day) == (2026, 1, 1)
    assert (stamped.hour, stamped.minute) == (8, 0)
    assert stamped.isoformat() == "2026-01-01T08:00:00+00:00"


def test_ingested_at_already_in_utc_is_stored_with_an_explicit_offset(store):
    """A UTC input keeps its offset on disk; naive-looking strings are corrupt."""
    when = datetime(2026, 7, 31, 12, 34, 56, tzinfo=timezone.utc)
    record_id = store.append(
        "s1", CALLS, guideline_version="cpic-2026-07", ingested_at=when
    )
    conn = raw(store)
    try:
        stored = conn.execute(
            "SELECT ingested_at FROM records WHERE record_id = ?", (record_id,)
        ).fetchone()[0]
    finally:
        conn.close()
    assert stored == "2026-07-31T12:34:56+00:00"


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


# --- Corrupt stored state is refused, not answered --------------------------
#
# Neither shape below is reachable through append(), so meeting one means the
# file was written by something that bypassed this module. The store's job then
# is to be loud, because both quiet answers are wrong in the dangerous
# direction: an empty call list reads as "no data for every gene", and a naive
# stamp reinterpreted as local time is a falsified audit trail.


def test_a_record_with_no_gene_calls_raises_rather_than_reading_as_no_data(store):
    """An empty result must not masquerade as "we looked and found nothing"."""
    record_id = insert_bare_record(store, "ghost")

    # The record really is there -- this is not the missing-record case.
    assert store.history("ghost") == [record_id]
    assert store.record_versions(record_id) == (PHARMCAT_VERSION, "cpic-2026-07")

    with pytest.raises(CorruptRecordError):
        store.record_calls(record_id)
    with pytest.raises(CorruptRecordError):
        store.latest("ghost")


def test_absent_and_corrupt_are_different_outcomes(store):
    """Three distinct outcomes, and they must not collapse into each other.

    * a subject with no records at all -> [] . A legitimate absence.
    * a subject whose newest record has no calls -> CorruptRecordError.
    * an unknown record_id -> KeyError.

    Collapsing the middle into the first is the failure this software must never
    make: presenting "the file is broken" as "this patient has no findings".
    """
    record_id = insert_bare_record(store, "ghost")

    assert store.latest("nobody-at-all") == []
    assert store.history("nobody-at-all") == []

    with pytest.raises(CorruptRecordError):
        store.latest("ghost")

    with pytest.raises(KeyError) as excinfo:
        store.record_calls(record_id + 12345)
    assert not isinstance(excinfo.value, CorruptRecordError)


def test_corrupt_record_error_is_a_value_error(store):
    """Callers already treating a bad record as bad input keep working."""
    record_id = insert_bare_record(store, "ghost")
    assert issubclass(CorruptRecordError, ValueError)
    with pytest.raises(ValueError):
        store.record_calls(record_id)


def test_a_corrupt_newest_record_does_not_hide_behind_an_older_good_one(store):
    """latest() must not silently fall back to the last readable record.

    Answering with the previous record would be worse than raising: the caller
    would get a real-looking phenotype attributed to the wrong record's
    versions, and would have no way to know.
    """
    good = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    corrupt = insert_bare_record(store, "s1")

    assert store.history("s1") == [good, corrupt]
    with pytest.raises(CorruptRecordError):
        store.latest("s1")
    # The older record is still readable on its own terms.
    assert store.record_calls(good) == sorted(CALLS, key=lambda c: c.gene)


def test_a_naive_stored_stamp_raises_instead_of_being_localized(store):
    """A stamp with no offset cannot be placed on a timeline without guessing.

    `datetime.fromisoformat("2020-01-01T00:00:00").astimezone(utc)` assumes the
    *reader's* zone, so the same immutable file reported a different ingest time
    depending on where it was opened -- on a machine in UTC-8 this stamp came
    back as 08:00 UTC. Raising is the only honest answer.
    """
    record_id = insert_bare_record(
        store, "s1", ingested_at="2020-01-01T00:00:00"
    )

    with pytest.raises(CorruptRecordError):
        store.record_ingested_at(record_id)


@pytest.mark.parametrize(
    "stamp",
    [
        pytest.param("2020-01-01T00:00:00", id="naive-datetime"),
        pytest.param("2020-01-01 00:00:00", id="naive-space-separated"),
        pytest.param("2020-01-01", id="naive-date-only"),
        pytest.param("not a timestamp at all", id="garbage"),
        pytest.param("2020-13-45T99:99:99+00:00", id="out-of-range-fields"),
        pytest.param("0", id="bare-zero"),
        pytest.param("2020-01-01T00:00:00+", id="truncated-offset"),
    ],
)
def test_an_unusable_stored_stamp_raises(store, stamp):
    record_id = insert_bare_record(store, "s1", ingested_at=stamp)
    with pytest.raises(CorruptRecordError):
        store.record_ingested_at(record_id)


def test_a_stamp_with_a_non_utc_offset_is_still_readable(store):
    """Only naive and unparseable stamps are corrupt; an offset is enough.

    A file written by an older version, or by another tool, that stored a
    correct instant in a different zone is not corrupt -- it is unambiguous. The
    reader normalizes it rather than rejecting it.
    """
    record_id = insert_bare_record(
        store, "s1", ingested_at="2026-01-01T00:00:00-08:00"
    )
    stamped = store.record_ingested_at(record_id)
    assert stamped == datetime(2026, 1, 1, 8, 0, tzinfo=timezone.utc)
    assert stamped.utcoffset() == timedelta(0)


def test_unknown_record_id_still_raises_key_error_not_corrupt_record(store):
    """A missing record and a broken record are different diagnoses."""
    for lookup in (
        lambda: store.record_calls(999),
        lambda: store.record_ingested_at(999),
        lambda: store.record_versions(999),
    ):
        with pytest.raises(KeyError) as excinfo:
            lookup()
        assert not isinstance(excinfo.value, CorruptRecordError)


# --- latest() is one consistent read ----------------------------------------


def test_latest_does_not_round_trip_through_history(store, monkeypatch):
    """latest() must resolve the record and its calls on a single connection.

    The race itself is not reachable from a single-threaded test: the old
    implementation called `history()` on one connection, then `record_calls()`
    on another, and a concurrent `append()` committing in between meant the
    calls returned could belong to a record newer than the id selected -- a
    result that never existed as a single state of the store.

    What *is* testable is the observable consequence of the fix: latest() no
    longer makes that second round-trip. Breaking `history()` therefore cannot
    break `latest()`. This is a proxy for the single-connection property, but a
    load-bearing one -- reverting latest() to `self.record_calls(self.history(
    subject_id)[-1])` fails here immediately.
    """
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    newest = [GeneCall("CYP2C19", "*2/*2", "Poor Metabolizer", CALLED)]
    second = store.append("s1", newest, guideline_version="cpic-2026-08")

    def exploded(self, subject_id):
        raise AssertionError(
            "latest() must not depend on history(): a second round-trip on a "
            "second connection is the race this fix removed"
        )

    monkeypatch.setattr(RecordStore, "history", exploded)

    assert store.latest("s1") == newest
    assert store.latest("s1")[0].phenotype == "Poor Metabolizer"
    assert store.latest("unknown-subject") == []

    # ...and the record it chose really is the most recently appended one.
    monkeypatch.undo()
    assert store.history("s1")[-1] == second


def test_latest_does_not_round_trip_through_record_calls(store, monkeypatch):
    """The other half: the calls are fetched on the connection already open."""
    store.append("s1", CALLS, guideline_version="cpic-2026-07")

    def exploded(self, record_id):
        raise AssertionError(
            "latest() must not reopen a connection via record_calls()"
        )

    monkeypatch.setattr(RecordStore, "record_calls", exploded)
    assert store.latest("s1") == sorted(CALLS, key=lambda c: c.gene)


def test_latest_opens_exactly_one_connection(store):
    """Pin the single-connection property directly by counting connects.

    monkeypatching history()/record_calls() shows latest() does not depend on
    those two methods; this shows the stronger thing they were a proxy for --
    that the whole read happens on one handle. Any reimplementation that splits
    the read across two connections, even without calling a public method,
    fails here.
    """
    store.append("s1", CALLS, guideline_version="cpic-2026-07")

    real_connect = sqlite3.connect
    opened = []

    def counting_connect(*args, **kwargs):
        opened.append(args[0] if args else kwargs.get("database"))
        return real_connect(*args, **kwargs)

    import pgxrecord.store as store_module

    original = store_module.sqlite3.connect
    store_module.sqlite3.connect = counting_connect
    try:
        calls = store.latest("s1")
    finally:
        store_module.sqlite3.connect = original

    assert calls == sorted(CALLS, key=lambda c: c.gene)
    assert len(opened) == 1, f"latest() opened {len(opened)} connections: {opened}"


def test_latest_returns_the_most_recently_appended_record_not_the_largest(store):
    """The contract in its own right, independent of how it is implemented."""
    ids = []
    for month, phenotype in enumerate(
        ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"], 1
    ):
        ids.append(
            store.append(
                "s1",
                [GeneCall("CYP2C19", f"*1/*{month}", phenotype, CALLED)],
                guideline_version=f"cpic-2026-{month:02d}",
            )
        )

    assert store.latest("s1") == store.record_calls(ids[-1])
    assert store.latest("s1")[0].phenotype == "Poor Metabolizer"
    # Another subject's later append must not become s1's latest.
    store.append("s2", CALLS, guideline_version="cpic-2026-09")
    assert store.latest("s1")[0].phenotype == "Poor Metabolizer"


def _open_handles(db_path):
    """How many of this process's fds point at `db_path`, by (device, inode).

    Counting fds is the only way to make the leak test able to fail. The
    previous version of this test just ran 59 appends and asserted the history
    length: deleting `finally: conn.close()` passed it, because 59 stray fds do
    not come close to the descriptor limit (1048576 here). Matching on the
    stat identity rather than on a readlink keeps this portable -- macOS
    /dev/fd/N is not readable as a symlink -- and survives the file being
    renamed or replaced mid-test.
    """
    target = os.stat(db_path)
    key = (target.st_dev, target.st_ino)
    count = 0
    for fd in range(4096):
        try:
            info = os.fstat(fd)
        except OSError:
            continue
        if (info.st_dev, info.st_ino) == key:
            count += 1
    return count


def test_open_handle_counter_can_actually_see_a_leak(tmp_path):
    """Guard the guard: prove _open_handles() responds to a real leak.

    Without this, a silently broken counter would make the leak test below
    unfalsifiable again -- the exact defect it was written to fix.
    """
    db_path = tmp_path / "records.db"
    RecordStore(db_path)
    baseline = _open_handles(db_path)

    leaked = [sqlite3.connect(db_path) for _ in range(4)]
    try:
        for conn in leaked:
            conn.execute("SELECT 1 FROM records LIMIT 1")
        assert _open_handles(db_path) == baseline + 4
    finally:
        for conn in leaked:
            conn.close()
    assert _open_handles(db_path) == baseline


def test_store_does_not_leak_connections_across_many_appends(tmp_path):
    """Every helper must close its handle; sqlite3's own context manager does not.

    `with sqlite3.connect(...)` commits but does not close, so a store built on
    it leaks one descriptor per call. Asserted by counting descriptors pointed
    at the database file, which is what makes this test capable of failing:
    removing `finally: conn.close()` from `_connect` leaves ~180 open handles
    here and the assertion catches it immediately.
    """
    db_path = tmp_path / "records.db"
    store = RecordStore(db_path)
    baseline = _open_handles(db_path)

    for month in range(1, 60):
        store.append("s1", CALLS, guideline_version=f"g-{month}")
        store.latest("s1")
        store.history("s1")
        store.record_versions(store.history("s1")[-1])

    assert len(store.history("s1")) == 59
    # No handle survives a completed call. Not "few", none: a store that leaks
    # one descriptor per read fails on a long-running process, and the whole
    # point of _connect's try/finally is that the leak is zero rather than slow.
    assert _open_handles(db_path) == baseline == 0
