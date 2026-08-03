"""Command line interface.

This is where the invariant becomes something a person reads, so the wording is
part of the contract rather than presentation:

* `cannot_assess` is printed as a loud CANNOT ASSESS banner and still carries
  the CPIC citation. It is never a quiet empty result and never borrows a phrase
  ("no interaction", "normal", "safe", "clear") that a reader scanning the
  output could take as an all-clear. See evaluate.py for why.
* A multi-gene answer also prints one OVERALL line, computed by
  `evaluate.overall_outcome`, which is the LEAST reassuring component. Warfarin
  has two pairs; reporting `guidance_found` because CYP2C9 was called while
  VKORC1 was never covered is the collapse this project exists to prevent,
  arriving one level above the code that forbids it.
* Every failure is surfaced. A pair table that will not load, a PharmCAT that
  will not run, a gene symbol the store rejects -- all of them exit non-zero
  with the message on stderr. The one thing this CLI must never do is print
  "nothing to report" because something upstream broke.

Exit codes are three-valued for the same reason the outcomes are:

    0  the query was answered (guidance_found or no_guidance_for_pair)
    1  the command failed and answered nothing
    2  the query could not be assessed -- an answer, but not a finding

2 exists so that a script wrapping this tool cannot treat "we do not know" as
success. `if pgx query codeine; then ...` reads as an all-clear on exit 0, so
cannot_assess must not exit 0.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from collections.abc import Sequence
from importlib.resources import files
from pathlib import Path

from pharmacogenomic_record import POSITIONS_FILENAME
from pharmacogenomic_record.caller import (
    GeneCall,
    PharmcatError,
    parse_phenotype_json,
    run_pharmcat,
)
from pharmacogenomic_record.drift import affected_by_guideline_change
from pharmacogenomic_record.evaluate import (
    CANNOT_ASSESS,
    GUIDANCE_FOUND,
    NO_GUIDANCE_FOR_PAIR,
    overall_outcome,
    query_drug,
)
from pharmacogenomic_record.guidelines import GuidelineTableError, load_pairs
from pharmacogenomic_record.ingest.raw import UnsupportedRawFile, parse_23andme
from pharmacogenomic_record.ingest.vcf import CoverageReport, build_vcf
from pharmacogenomic_record.positions import load_positions
from pharmacogenomic_record.store import CorruptRecordError, RecordStore

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_CANNOT_ASSESS = 2

# The reference tables ship inside the package and are located through
# importlib.resources, not a filesystem walk relative to the source tree. The
# tool's answer must not depend on where the shell happens to be -- but it also
# must not depend on how the package was installed. A `parents[2] / "data"` walk
# only finds the tables in an editable checkout; a wheel installed into
# site-packages has no repo root two levels up, so the console script silently
# broke for every non-editable install. `files()` resolves the same package
# data whether the code runs from a checkout or an installed wheel.
_DATA_DIR = files("pharmacogenomic_record") / "data"
PAIRS_PATH = Path(str(_DATA_DIR / "gene_drug_pairs.json"))
POSITIONS_PATH = Path(str(_DATA_DIR / POSITIONS_FILENAME))

# The version stamp written onto every record, naming the revision of the pair
# table that produced it. Pinned here beside the table it describes: a record
# that does not know which guidance produced it cannot be re-evaluated when that
# guidance moves, which is the entire purpose of storing it.
GUIDELINE_VERSION = "cpic-2026-07"

# How each outcome is labelled. Spelled once so the banner a user greps for and
# the OVERALL summary cannot drift apart.
_LABELS = {
    CANNOT_ASSESS: "CANNOT ASSESS",
    GUIDANCE_FOUND: "GUIDANCE",
    NO_GUIDANCE_FOR_PAIR: "NO CPIC PAIR",
}

# Printed under every query. Deliberately worded without any of the phrases the
# reassurance checks forbid: a disclaimer that reads as comfort is worse than
# none. "Reference" and "citation" are the claims this tool can support.
_DISCLAIMER = (
    "This output is a reference to published CPIC citations, keyed on stored "
    "gene calls. It is not a medical device, not clinical decision support, and "
    "not a basis for any treatment decision."
)


def ingest_to_calls(
    raw_path: Path, positions_path: Path, workdir: Path
) -> tuple[Path, CoverageReport]:
    """Parse a raw file and write a VCF. Does not require Docker.

    Split out from `cmd_ingest` because it is the half that can be tested: it
    has no Docker dependency, and it is where every coverage decision is made.
    """
    calls = parse_23andme(raw_path)
    positions = load_positions(positions_path)
    workdir.mkdir(parents=True, exist_ok=True)
    vcf_path = workdir / f"{raw_path.stem}.vcf"
    report = build_vcf(calls, positions, vcf_path)
    return vcf_path, report


def calls_from_phenotype(
    phenotype_json: Path, report: CoverageReport
) -> list[GeneCall]:
    """Translate PharmCAT output under the coverage report's strict rule.

    The only place the two modules meet, and the reason it is a named function
    rather than two keyword arguments at the call site: the three sets have to
    stay coherent, and `parse_phenotype_json` cannot check that for us.

    `genes_fully_covered` is passed by OMISSION -- it is the only set whose
    members are allowed to reach PharmCAT's own answer. Passing
    `genes_partially_covered` where the older, looser field of the same name
    once lived would send every gene with a single covered position through as
    eligible to be called, so the mapping is stated explicitly here:

        genes_fully_uncovered    -> uncovered_genes        -> not_covered
        genes_partially_covered  -> partially_covered_genes -> indeterminate
        genes_fully_covered      -> neither                -> PharmCAT decides
    """
    return parse_phenotype_json(
        phenotype_json,
        uncovered_genes=report.genes_fully_uncovered,
        partially_covered_genes=report.genes_partially_covered,
    )


def cmd_ingest(
    store: RecordStore,
    raw_path: Path,
    subject_id: str,
    workdir: Path,
    positions_path: Path = POSITIONS_PATH,
) -> None:
    """Ingest a raw export: VCF, PharmCAT, store.

    Coverage is printed before PharmCAT is invoked, so that a run which fails at
    the Docker step has still told the user what their array actually covered.
    Nothing is caught here: a PharmCAT failure must reach `main` and exit
    non-zero rather than leave a partial record behind.
    """
    vcf_path, report = ingest_to_calls(raw_path, positions_path, workdir)
    print(f"wrote {vcf_path}")
    print(f"covered positions: {len(report.covered_rsids)}")
    print(f"uncovered positions: {len(report.uncovered_rsids)}")
    print(
        f"positions with no rsID, unjoinable from a 23andMe file: "
        f"{report.unjoinable_positions}"
    )
    print(
        f"genes fully covered (eligible to be called): "
        f"{len(report.genes_fully_covered)}"
        + (
            f" -- {', '.join(sorted(report.genes_fully_covered))}"
            if report.genes_fully_covered
            else ""
        )
    )
    print(
        f"genes partially covered (recorded indeterminate): "
        f"{len(report.genes_partially_covered)}"
        + (
            f" -- {', '.join(sorted(report.genes_partially_covered))}"
            if report.genes_partially_covered
            else ""
        )
    )
    print(
        f"genes with no coverage (recorded not_covered): "
        f"{len(report.genes_fully_uncovered)}"
        + (
            f" -- {', '.join(sorted(report.genes_fully_uncovered))}"
            if report.genes_fully_uncovered
            else ""
        )
    )

    phenotype_json = run_pharmcat(vcf_path, workdir)
    gene_calls = calls_from_phenotype(phenotype_json, report)
    record_id = store.append(
        subject_id, gene_calls, guideline_version=GUIDELINE_VERSION
    )
    print(f"stored record {record_id} for subject {subject_id}")


def cmd_query(
    store: RecordStore,
    subject_id: str,
    drug: str,
    pairs_path: Path = PAIRS_PATH,
) -> str:
    """Print one line per relevant gene-drug pair, then an overall verdict.

    Returns the overall outcome so `main` can choose an exit code without
    re-deriving it. The per-pair lines are printed in full even when the summary
    is CANNOT ASSESS: the summary tells the user how much to trust the answer,
    it does not replace the answer.
    """
    pairs = load_pairs(pairs_path)
    results = query_drug(store, subject_id, drug, pairs)

    for result in results:
        print(f"[{_LABELS[result.outcome]}] {result.explanation}")

    overall = overall_outcome(results)
    print(f"OVERALL: {_LABELS[overall]}")
    print(_DISCLAIMER)
    return overall


def cmd_drift(
    store: RecordStore,
    changed_pair_ids: Sequence[str],
    pairs_path: Path = PAIRS_PATH,
) -> None:
    """Report which stored records a guideline revision touches.

    `changed_pair_ids` is a Sequence, and `argparse` builds it with
    `action="append"`, so a single `--changed-pair` still arrives as a
    one-element list. `affected_by_guideline_change` raises on a bare string on
    purpose -- iterating one would compare its characters and report that the
    revision affected nobody.
    """
    pairs = load_pairs(pairs_path)
    affected = affected_by_guideline_change(store, list(changed_pair_ids), pairs)

    ids = ", ".join(changed_pair_ids)
    if not affected:
        # An explicit negative, and scoped to what it actually establishes. This
        # says which records the revision touches; it says nothing about whether
        # the revision matters to anyone whose genotype we have never seen.
        print(
            f"No stored record holds a gene belonging to the changed pair(s) "
            f"{ids}. That states only which stored records this revision "
            f"touches."
        )
    else:
        for record in affected:
            print(
                f"[AFFECTED] subject {record.subject_id} gene {record.gene} "
                f"pair(s) {', '.join(record.changed_pair_ids)}"
            )
    print(_DISCLAIMER)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pharmacogenomic-record",
        description=(
            "Longitudinal pharmacogenomic record over PharmCAT. Reference "
            "tooling only: not a medical device."
        ),
    )
    parser.add_argument("--db", type=Path, default=Path("records.db"))
    sub = parser.add_subparsers(dest="command", required=True)

    p_ingest = sub.add_parser("ingest", help="ingest a 23andMe raw file")
    p_ingest.add_argument("raw_path", type=Path)
    p_ingest.add_argument("--subject", required=True)
    p_ingest.add_argument("--workdir", type=Path, default=Path("work"))
    p_ingest.add_argument("--positions", type=Path, default=POSITIONS_PATH)

    p_query = sub.add_parser("query", help="query a drug against a stored record")
    p_query.add_argument("drug")
    p_query.add_argument("--subject", required=True)
    p_query.add_argument("--pairs", type=Path, default=PAIRS_PATH)

    p_drift = sub.add_parser(
        "drift", help="report stored records a guideline revision touches"
    )
    p_drift.add_argument(
        "--changed-pair",
        action="append",
        required=True,
        dest="changed_pairs",
        metavar="CPIC_PAIR_ID",
        help="repeatable; a cpic_pair_id whose guidance changed",
    )
    p_drift.add_argument("--pairs", type=Path, default=PAIRS_PATH)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run one subcommand, returning its exit code.

    Every expected failure is caught here and only here, so that each one exits
    EXIT_ERROR with its message on stderr. The catch list is explicit rather
    than a bare `except Exception`: an unanticipated exception should still
    reach the user as a traceback, because a message this code did not write is
    not a message this code can promise is honest.

    `ValueError` is on the list because it is the vocabulary this codebase
    refuses in -- an empty pair table, a blank drug name, a naive timestamp --
    and `CorruptRecordError` and `GuidelineTableError` are both subclasses of it.
    The two are still named explicitly, for the reader rather than the
    interpreter. `TypeError` is not caught: nothing reachable through argparse
    raises one (drift's bare-string guard cannot fire, since `action="append"`
    always yields a list), so a TypeError here means a bug, and a bug deserves a
    traceback rather than a tidy one-line "error:".
    """
    args = _build_parser().parse_args(argv)

    try:
        store = RecordStore(args.db)
        if args.command == "ingest":
            cmd_ingest(
                store,
                args.raw_path,
                args.subject,
                args.workdir,
                positions_path=args.positions,
            )
            return EXIT_OK
        if args.command == "drift":
            cmd_drift(store, args.changed_pairs, pairs_path=args.pairs)
            return EXIT_OK

        outcome = cmd_query(store, args.subject, args.drug, pairs_path=args.pairs)
        # guidance_found and no_guidance_for_pair are both answers about the
        # table and the record. cannot_assess is not, and must not be reported
        # to a calling script as success.
        return EXIT_CANNOT_ASSESS if outcome == CANNOT_ASSESS else EXIT_OK
    except (
        UnsupportedRawFile,
        PharmcatError,
        GuidelineTableError,
        CorruptRecordError,
        sqlite3.Error,
        OSError,
        KeyError,
        ValueError,
    ) as err:
        # stdout is flushed first, deliberately. When stdout is redirected to a
        # file or a pipe it is block-buffered while stderr is not, so without
        # this the error line lands *above* output that was produced before it --
        # e.g. "cannot run PharmCAT" printed above the coverage summary, making
        # the failure look like it happened earlier than it did. On a terminal
        # the order is already right; this makes it right everywhere.
        sys.stdout.flush()
        print(f"error: {err}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":  # pragma: no cover - exercised via subprocess
    raise SystemExit(main())
