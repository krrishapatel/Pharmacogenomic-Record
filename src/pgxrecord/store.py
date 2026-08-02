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

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import CALLED, INDETERMINATE, NOT_COVERED, GeneCall

# The closed coverage vocabulary, taken from pgxrecord.caller rather than
# re-spelled here: two copies of a vocabulary is one copy too many.
_COVERAGE_VALUES = (CALLED, NOT_COVERED, INDETERMINATE)
_COVERAGE_LIST = ", ".join(f"'{value}'" for value in _COVERAGE_VALUES)
_UNCALLED_LIST = ", ".join(f"'{value}'" for value in (NOT_COVERED, INDETERMINATE))

_SCHEMA = f"""
CREATE TABLE IF NOT EXISTS records (
    record_id INTEGER PRIMARY KEY AUTOINCREMENT,
    subject_id TEXT NOT NULL CHECK (trim(subject_id) <> ''),
    -- A record that does not know which tool and guideline produced it is
    -- worthless, so the stamp is NOT NULL *and* non-blank: sqlite treats ''
    -- as a perfectly good NOT NULL value.
    pharmcat_version TEXT NOT NULL CHECK (trim(pharmcat_version) <> ''),
    guideline_version TEXT NOT NULL CHECK (trim(guideline_version) <> ''),
    ingested_at TEXT NOT NULL CHECK (trim(ingested_at) <> '')
);

CREATE TABLE IF NOT EXISTS gene_calls (
    record_id INTEGER NOT NULL REFERENCES records(record_id),
    gene TEXT NOT NULL CHECK (trim(gene) <> ''),
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
        OR (diplotype IS NOT NULL AND trim(diplotype) <> '')
    ),
    PRIMARY KEY (record_id, gene)
);

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
        """Gene calls from the subject's most recently appended record."""
        ids = self.history(subject_id)
        if not ids:
            return []
        return self.record_calls(ids[-1])

    def record_calls(self, record_id: int) -> list[GeneCall]:
        """The gene calls stored in one record, by gene name.

        Every field comes straight out of the row, so a stored
        not_covered/indeterminate row cannot be reassembled into a call: there
        is no default and no fallback to substitute for a NULL diplotype.
        """
        with self._connect() as conn:
            self._require_record(conn, record_id)
            rows = conn.execute(
                "SELECT gene, diplotype, phenotype, coverage FROM gene_calls "
                "WHERE record_id = ? ORDER BY gene",
                (record_id,),
            ).fetchall()
        return [GeneCall(*row) for row in rows]

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
        return datetime.fromisoformat(row[0]).astimezone(timezone.utc)

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
    def _require_record(conn: sqlite3.Connection, record_id: int) -> None:
        """Distinguish "no such record" from "a record with no calls".

        The latter cannot exist -- append() rejects an empty call list -- so an
        empty result means the record_id is wrong, and that is an error rather
        than an absence of findings.
        """
        row = conn.execute(
            "SELECT 1 FROM records WHERE record_id = ?", (record_id,)
        ).fetchone()
        if row is None:
            raise KeyError(f"no record {record_id}")
