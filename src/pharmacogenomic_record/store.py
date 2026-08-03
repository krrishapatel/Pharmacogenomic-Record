"""Append-only pharmacogenomic record store.

Immutability is a correctness requirement, not a preference. A phenotype call
only means something relative to the tool and guideline versions that produced
it. Overwriting a call destroys the ability to explain why the system once
said something different -- which is exactly the question that matters when
guidance is revised.

SQLite triggers enforce this so that a future careless UPDATE fails loudly
rather than silently rewriting clinical history. The application never issues
an UPDATE or a DELETE, but the guarantee does not rest on that: a future
caller who opens the file directly is stopped by the schema.

The other half of the guarantee is that nothing may be stored in a state that
reads back as more certain than it was written. That is why the coverage
vocabulary is a CHECK constraint rather than a convention, and why a
not_covered or indeterminate row is structurally forbidden a diplotype or a
phenotype: sqlite's dynamic typing will happily accept a value of the wrong
type in any column, so the constraints -- not the declared types -- are what
keep "we do not know" from decaying into "normal metabolizer".
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from pharmacogenomic_record import PHARMCAT_VERSION
from pharmacogenomic_record.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall

# The closed coverage vocabulary, taken from pharmacogenomic_record.caller rather than
# re-spelled here: two copies of a vocabulary is one copy too many.
_COVERAGE_VALUES = (CALLED, NOT_COVERED, INDETERMINATE)
_COVERAGE_LIST = ", ".join(f"'{value}'" for value in _COVERAGE_VALUES)
_UNCALLED_LIST = ", ".join(f"'{value}'" for value in (NOT_COVERED, INDETERMINATE))

# sqlite's one-argument trim() strips SPACES ONLY: `trim(char(9)) <> ''` is 1,
# so a `trim(x) <> ''` CHECK happily accepts a tab, a newline or a vertical tab
# as a subject id or a guideline version. A record stamped with a tab as its
# guideline version does not know which guidance produced it, which is exactly
# what the constraint exists to prevent. Every blank check therefore passes an
# explicit character set, spelled once here so the six uses cannot drift apart.
_BLANK_CHARS = "' ' || char(9) || char(10) || char(13) || char(11) || char(12)"


def _non_blank(column: str) -> str:
    """SQL that is true only when `column` holds a non-whitespace character."""
    return f"trim({column}, {_BLANK_CHARS}) <> ''"


_SCHEMA = f"""
CREATE TABLE IF NOT EXISTS records (
    record_id INTEGER PRIMARY KEY AUTOINCREMENT,
    subject_id TEXT NOT NULL CHECK ({_non_blank("subject_id")}),
    -- A record that does not know which tool and guideline produced it is
    -- worthless, so the stamp is NOT NULL *and* non-blank: sqlite treats ''
    -- as a perfectly good NOT NULL value.
    pharmcat_version TEXT NOT NULL CHECK ({_non_blank("pharmcat_version")}),
    guideline_version TEXT NOT NULL CHECK ({_non_blank("guideline_version")}),
    ingested_at TEXT NOT NULL CHECK ({_non_blank("ingested_at")})
);

-- WITHOUT ROWID is load-bearing, not a storage optimization. In an ordinary
-- rowid table, `INSERT OR REPLACE` that names an explicit rowid resolves the
-- conflict through an internal delete that skips BEFORE DELETE (sqlite only
-- fires it under `PRAGMA recursive_triggers`, which is per-connection and off
-- by default), and because the replacement row can carry a *different*
-- (record_id, gene) pair, the no-replace trigger below never fires either. The
-- net effect was that a stored phenotype could be deleted and a fabricated
-- gene put in its place, past all six triggers. Removing the implicit rowid
-- removes the only key that could be targeted that way: `records` is immune
-- for the same reason from the other direction -- record_id *is* its rowid, so
-- its primary-key trigger sees the conflict.
CREATE TABLE IF NOT EXISTS gene_calls (
    record_id INTEGER NOT NULL REFERENCES records(record_id),
    -- Gene symbols are joined across layers by exact string equality:
    -- `subjects_with_gene` compares with `g.gene = ?` and `query_drug` looks the
    -- gene up as a dict key, while the guideline table is force-uppercased by
    -- guidelines.normalize_gene. A stored "cyp2c19" or " CYP2C19" therefore
    -- matches nothing: the subject is not reported as affected, drops out of
    -- drift reports, and reads as "no data for this gene" -- a confident wrong
    -- answer with no error anywhere. HGNC symbols are uppercase by definition,
    -- so the canonical form is checkable rather than a guess.
    --
    -- Rejected here rather than upcased in Python on the way in, deliberately.
    -- Every PX= tag in the pinned PharmCAT positions table is uppercase; a
    -- lowercase symbol arriving means that property of an external tool's
    -- output has changed, and silently folding it would hide the version skew
    -- while this code kept resting on the assumption. Failing at the write is
    -- the only place the truth is still recoverable.
    --
    -- Two caveats this expression lives with, neither of which weakens it for
    -- real symbols (ASCII alphanumerics plus '-'):
    --   * sqlite's upper() folds ASCII only, so a non-ASCII symbol satisfies
    --     `gene = upper(gene)` yet still fails to match Python's str.upper()
    --     downstream. The CHECK narrows the hole, it does not close it.
    --   * a CHECK that evaluates to NULL counts as satisfied, and
    --     `NULL = upper(NULL)` is NULL -- so this clause alone would admit a
    --     NULL gene. NOT NULL above is what actually rejects it, and is load
    --     bearing for that reason, not decoration.
    -- Leading/trailing whitespace is pinned too: " CYP2C19" is invisible to the
    -- same lookups for the same reason, and multi-character trim() is used
    -- because the one-argument form strips spaces only (see _BLANK_CHARS).
    -- Internal whitespace is left alone -- no real symbol has any, and
    -- rejecting a legitimate symbol would be worse than the bug.
    gene TEXT NOT NULL CHECK (
        {_non_blank("gene")}
        AND gene = upper(gene)
        AND gene = trim(gene, {_BLANK_CHARS})
    ),
    diplotype TEXT,
    phenotype TEXT,
    coverage TEXT NOT NULL CHECK (coverage IN ({_COVERAGE_LIST})),
    -- A gene we did not call may not carry a call. Without this, a stray
    -- INSERT could store coverage='not_covered' alongside a diplotype and a
    -- phenotype, and any reader that looks at the phenotype column would
    -- report a confident result for a gene the array never informed.
    CHECK (
        coverage NOT IN ({_UNCALLED_LIST})
        OR (diplotype IS NULL AND phenotype IS NULL)
    ),
    -- Conversely, a "called" gene with nothing to report is a contradiction.
    -- The phenotype may legitimately be NULL (F2, F5, VKORC1 and friends have
    -- no metabolizer phenotype at all), but the diplotype may not.
    --
    -- The `IS NOT NULL` is load-bearing and not redundant with the trim(): a
    -- CHECK passes when it evaluates to NULL, and `trim(NULL) <> ''` is NULL,
    -- so `coverage <> 'called' OR trim(diplotype) <> ''` would silently admit
    -- exactly the row it is meant to reject -- coverage='called' with no
    -- diplotype at all. SQL three-valued logic defaults to permitting.
    CHECK (
        coverage <> '{CALLED}'
        OR (diplotype IS NOT NULL AND {_non_blank("diplotype")})
    ),
    PRIMARY KEY (record_id, gene)
) WITHOUT ROWID;

CREATE INDEX IF NOT EXISTS idx_gene_calls_gene ON gene_calls(gene);
CREATE INDEX IF NOT EXISTS idx_records_subject ON records(subject_id);

CREATE TRIGGER IF NOT EXISTS gene_calls_no_update
BEFORE UPDATE ON gene_calls
BEGIN
    SELECT RAISE(ABORT, 'gene_calls is append-only');
END;

CREATE TRIGGER IF NOT EXISTS gene_calls_no_delete
BEFORE DELETE ON gene_calls
BEGIN
    SELECT RAISE(ABORT, 'gene_calls is append-only');
END;

CREATE TRIGGER IF NOT EXISTS records_no_update
BEFORE UPDATE ON records
BEGIN
    SELECT RAISE(ABORT, 'records is append-only');
END;

CREATE TRIGGER IF NOT EXISTS records_no_delete
BEFORE DELETE ON records
BEGIN
    SELECT RAISE(ABORT, 'records is append-only');
END;

-- INSERT OR REPLACE is an overwrite wearing an INSERT's clothes: sqlite
-- implements it as delete-then-insert, but it SKIPS the BEFORE DELETE trigger
-- unless `PRAGMA recursive_triggers` is on -- and that pragma is
-- per-connection and off by default, so an outside caller gets the unsafe
-- setting for free. Without these BEFORE INSERT triggers, a single REPLACE
-- statement silently rewrites a stored phenotype past all four triggers
-- above. Rejecting an insert onto an existing key closes that door from the
-- schema side, where the guarantee belongs.
CREATE TRIGGER IF NOT EXISTS gene_calls_no_replace
BEFORE INSERT ON gene_calls
WHEN EXISTS (
    SELECT 1 FROM gene_calls
    WHERE record_id = NEW.record_id AND gene = NEW.gene
)
BEGIN
    SELECT RAISE(ABORT, 'gene_calls is append-only');
END;

CREATE TRIGGER IF NOT EXISTS records_no_replace
BEFORE INSERT ON records
WHEN EXISTS (SELECT 1 FROM records WHERE record_id = NEW.record_id)
BEGIN
    SELECT RAISE(ABORT, 'records is append-only');
END;
"""


class CorruptRecordError(ValueError):
    """A stored record cannot be read back as what it claims to be.

    Raised rather than returning a plausible-looking answer. The two shapes it
    covers -- a record with no gene calls, and a stamp that is not an aware
    timestamp -- are both unreachable through `append()`, so encountering one
    means the file was written by something that bypassed this module. In both
    cases the honest answer is louder than the convenient one: an empty call
    list is indistinguishable from "no data for every gene", and a naive stamp
    silently reinterpreted as local time is a falsified audit trail.

    Subclasses ValueError so that callers already treating a bad record as bad
    input keep working.
    """


class RecordStore:
    """Append-only store of pharmacogenomic records."""

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        with self._connect() as conn:
            conn.executescript(_SCHEMA)

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        """Open a connection, commit on success, roll back on error, always close.

        sqlite3's own connection context manager commits the transaction but
        does NOT close the handle, so using `with sqlite3.connect(...)`
        directly leaks a file descriptor per call. This wraps both: the inner
        `with conn` gives commit-or-rollback, the outer try/finally gives close.
        """
        conn = sqlite3.connect(self.db_path)
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            with conn:
                yield conn
        finally:
            conn.close()

    def append(
        self,
        subject_id: str,
        calls: Sequence[GeneCall],
        guideline_version: str,
        ingested_at: datetime | None = None,
    ) -> int:
        """Write a new record. Never modifies an existing one.

        Re-ingesting a subject appends; it does not replace. The whole write is
        one transaction, so a failure anywhere -- a bad coverage value, a
        duplicated gene -- leaves no record at all rather than a record whose
        gene list is a prefix of the truth.

        `ingested_at` must be timezone-aware and is normalized to UTC. It is a
        parameter rather than a buried `datetime.now()` so that callers
        replaying an archive, and tests asserting on the stamp, can both say
        exactly when the observation happened.
        """
        if not calls:
            raise ValueError(
                f"refusing to store an empty record for {subject_id!r}: a record "
                f"with no gene calls reads back as 'no data' for every gene"
            )

        if ingested_at is None:
            ingested_at = datetime.now(timezone.utc)
        if ingested_at.tzinfo is None or ingested_at.utcoffset() is None:
            raise ValueError(
                f"ingested_at must be timezone-aware, got {ingested_at!r}; a "
                f"naive timestamp cannot be compared across ingests"
            )
        stamp = ingested_at.astimezone(timezone.utc).isoformat()

        with self._connect() as conn:
            cursor = conn.execute(
                "INSERT INTO records "
                "(subject_id, pharmcat_version, guideline_version, ingested_at) "
                "VALUES (?, ?, ?, ?)",
                (subject_id, PHARMCAT_VERSION, guideline_version, stamp),
            )
            record_id = int(cursor.lastrowid)
            conn.executemany(
                "INSERT INTO gene_calls "
                "(record_id, gene, diplotype, phenotype, coverage) "
                "VALUES (?, ?, ?, ?, ?)",
                [
                    (record_id, c.gene, c.diplotype, c.phenotype, c.coverage)
                    for c in calls
                ],
            )
        return record_id

    def history(self, subject_id: str) -> list[int]:
        """Every record_id for a subject, oldest first.

        Ordered by record_id, i.e. by the order the records were appended.
        That is deliberately insertion order and not `ingested_at` order: a
        backdated replay must not be able to rewrite what the store considers
        its most recent knowledge.
        """
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT record_id FROM records WHERE subject_id = ? "
                "ORDER BY record_id",
                (subject_id,),
            ).fetchall()
        return [row[0] for row in rows]

    def latest(self, subject_id: str) -> list[GeneCall]:
        """Gene calls from the subject's most recently appended record.

        One connection for both halves of the read. Choosing the record on one
        connection and then fetching its calls on another leaves a window in
        which a concurrent `append()` commits, and the calls that come back can
        belong to a record newer than the id this call selected -- i.e. a result
        that never existed as a single consistent state of the store.
        """
        with self._connect() as conn:
            row = conn.execute(
                "SELECT record_id FROM records WHERE subject_id = ? "
                "ORDER BY record_id DESC LIMIT 1",
                (subject_id,),
            ).fetchone()
            if row is None:
                return []
            return self._calls_for_record(conn, int(row[0]))

    def record_calls(self, record_id: int) -> list[GeneCall]:
        """The gene calls stored in one record, by gene name.

        Every field comes straight out of the row, so a stored
        not_covered/indeterminate row cannot be reassembled into a call: there
        is no default and no fallback to substitute for a NULL diplotype.
        """
        with self._connect() as conn:
            return self._calls_for_record(conn, record_id)

    def record_versions(self, record_id: int) -> tuple[str, str]:
        """The (pharmcat_version, guideline_version) a record was made with."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT pharmcat_version, guideline_version FROM records "
                "WHERE record_id = ?",
                (record_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"no record {record_id}")
        return (row[0], row[1])

    def record_ingested_at(self, record_id: int) -> datetime:
        """When a record was ingested, as a UTC-aware datetime."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT ingested_at FROM records WHERE record_id = ?",
                (record_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"no record {record_id}")

        stamp = row[0]
        try:
            parsed = datetime.fromisoformat(stamp)
        except (TypeError, ValueError) as err:
            raise CorruptRecordError(
                f"record {record_id} has an unparseable ingested_at {stamp!r}: "
                f"{err}"
            ) from err
        # A naive stored stamp must not be quietly localized. astimezone() on a
        # naive datetime assumes the *reader's* local zone, so the same file
        # would report a different ingest time depending on where it was
        # opened -- a falsified audit trail on a record that is supposed to be
        # immutable. append() only ever writes an offset, so a naive stamp means
        # something bypassed this module.
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            raise CorruptRecordError(
                f"record {record_id} has a naive ingested_at {stamp!r}; a stored "
                f"stamp without an offset cannot be placed on a timeline "
                f"without guessing a zone"
            )
        return parsed.astimezone(timezone.utc)

    def subjects_with_gene(self, gene: str) -> list[str]:
        """Subjects with any stored call for this gene, across all records.

        Deliberately not restricted to the latest record: a guideline change
        is relevant to anyone ever genotyped for the gene.
        """
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT DISTINCT r.subject_id FROM records r "
                "JOIN gene_calls g ON g.record_id = r.record_id "
                "WHERE g.gene = ? ORDER BY r.subject_id",
                (gene,),
            ).fetchall()
        return [row[0] for row in rows]

    @staticmethod
    def _calls_for_record(
        conn: sqlite3.Connection, record_id: int
    ) -> list[GeneCall]:
        """Read one record's calls, refusing to return an empty list.

        Three outcomes, kept distinct on purpose:

        * no `records` row -> KeyError. The record_id is wrong; that is a bug in
          the caller, not an absence of findings.
        * a `records` row with no `gene_calls` -> CorruptRecordError. An empty
          list would be read downstream as "no data for every gene", which is
          the collapse the coverage vocabulary exists to prevent. `append()`
          rejects an empty call list, so this shape can only arrive via raw SQL
          -- and when it does, the store must say so rather than answer.
        * a `records` row with calls -> the calls.
        """
        row = conn.execute(
            "SELECT 1 FROM records WHERE record_id = ?", (record_id,)
        ).fetchone()
        if row is None:
            raise KeyError(f"no record {record_id}")

        rows = conn.execute(
            "SELECT gene, diplotype, phenotype, coverage FROM gene_calls "
            "WHERE record_id = ? ORDER BY gene",
            (record_id,),
        ).fetchall()
        if not rows:
            raise CorruptRecordError(
                f"record {record_id} is corrupt: the record exists but stores no "
                f"gene calls, and an empty result would read as 'no data' for "
                f"every gene rather than as a missing record"
            )
        return [GeneCall(*call) for call in rows]
