"""CLI tests.

Two things are pinned here that exist nowhere else:

1. The strict coverage rule. A gene is eligible to be `called` only when EVERY
   rsID-joinable position for that gene was covered by the array. One position
   short makes it indeterminate. Loosening that to "at least one position
   covered" is what would let PharmCAT's reference-assumption at unobserved
   positions be reported as a confident diplotype.
2. The rendered wording of the three outcomes. `cannot_assess` must arrive as a
   banner that cites the guideline and contains no phrase a reader (or a
   string-matching consumer) could take as an all-clear.

Every path is anchored on `Path(__file__)`, never on the CWD: the suite is run
from the repo root and from /tmp and must behave identically.
"""

import json
import sys
from pathlib import Path

import pytest

from pharmacogenomic_record import POSITIONS_FILENAME
from pharmacogenomic_record.caller import (
    CALLED,
    INDETERMINATE,
    NOT_COVERED,
    GeneCall,
    PharmcatError,
)
from pharmacogenomic_record.cli import (
    EXIT_CANNOT_ASSESS,
    EXIT_ERROR,
    EXIT_OK,
    calls_from_phenotype,
    cmd_drift,
    cmd_ingest,
    cmd_query,
    ingest_to_calls,
    main,
)
from pharmacogenomic_record.evaluate import CANNOT_ASSESS, GUIDANCE_FOUND
from pharmacogenomic_record.ingest.raw import RawCall
from pharmacogenomic_record.ingest.vcf import build_vcf
from pharmacogenomic_record.positions import ReferencePosition
from pharmacogenomic_record.store import RecordStore

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
DATA = REPO_ROOT / "src" / "pharmacogenomic_record" / "data"
POSITIONS = DATA / POSITIONS_FILENAME
PAIRS = DATA / "gene_drug_pairs.json"
PHENOTYPE_SAMPLE = FIXTURES / "pharmcat_phenotype_sample.json"

# Phrases that would turn "we do not know" into an all-clear. Same list Task 7
# applies to the explanations; applied here to the *rendered* output, because a
# safe explanation printed under a reassuring banner is still reassuring.
REASSURING = (
    "no interaction",
    "no issue",
    "safe",
    "normal",
    "no guidance",
    "clear",
    "fine",
    "nothing found",
    "no data for",
)

# A subject id is caller-supplied text, so it can contain every phrase above.
# It must never be echoed into the answer -- see evaluate.query_drug.
HOSTILE_SUBJECT_ID = "patient A - no interaction, safe, normal"

# Two joinable CYP2C19 positions and one joinable DPYD position, so a gene can
# be taken from "one position short" to "every position covered" by adding a
# single raw call. Coordinates are real 3.4.0 values; the point being tested is
# the counting rule, not the coordinates.
#
# The final entry is synthetic: a gene reachable ONLY through an rsID-less
# position, which is the vacuous case of the strict rule. It models the join
# mechanism and is not a claim about real CFTR data.
REF = [
    ReferencePosition(
        chrom="chr10", pos=94761900, rsid="rs12248560", ref="C", alt=("T",),
        gene="CYP2C19",
    ),
    ReferencePosition(
        chrom="chr10", pos=94781859, rsid="rs4244285", ref="G", alt=("A",),
        gene="CYP2C19",
    ),
    ReferencePosition(
        chrom="chr1", pos=97544543, rsid="rs114096998", ref="G", alt=("T",),
        gene="DPYD",
    ),
    ReferencePosition(
        chrom="chr7", pos=117509089, rsid=None, ref="A", alt=("G",), gene="CFTR"
    ),
]

CYP2C19_CALLS = [
    RawCall(rsid="rs12248560", chrom="10", pos=1, genotype="CC"),
    RawCall(rsid="rs4244285", chrom="10", pos=2, genotype="AG"),
]


def assert_no_reassurance(text: str) -> None:
    lowered = text.lower()
    for phrase in REASSURING:
        assert phrase not in lowered, f"output reads as an all-clear: {phrase!r}"


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


# --------------------------------------------------------------------------
# The reference tables are located as package data, not by a filesystem walk.
# --------------------------------------------------------------------------


def test_default_data_paths_resolve_through_the_package_not_a_source_walk():
    """The pinned tables must be found via the package, not `parents[2]/data`.

    A `Path(__file__).resolve().parents[2] / "data"` walk finds the tables only
    in an editable checkout: a wheel installed into site-packages has no repo
    root two levels up, so `query` failed to load its pair table for every
    non-editable install. Resolving through `importlib.resources` finds the
    same data whether the code runs from a checkout or an installed wheel.

    Anchored on the package location, so it fails if the defaults ever revert
    to walking up from the source file.
    """
    import pharmacogenomic_record as pkg
    from pharmacogenomic_record import cli

    package_data = Path(pkg.__file__).resolve().parent / "data"

    assert cli.PAIRS_PATH == package_data / "gene_drug_pairs.json"
    assert cli.POSITIONS_PATH == package_data / POSITIONS_FILENAME
    assert cli.PAIRS_PATH.is_file()
    assert cli.POSITIONS_PATH.is_file()
    # Readable through the resolved path, which is the property `query` needs.
    assert cli.PAIRS_PATH.read_text(encoding="utf-8").strip()
    assert cli.POSITIONS_PATH.read_text(encoding="utf-8").strip()


def test_default_data_paths_are_independent_of_the_working_directory(tmp_path):
    """Resolved from an unrelated CWD, the defaults must still find the tables.

    Run in a subprocess whose CWD is a scratch directory with no `data/` in
    sight, so a resolver that depended on the working directory -- or on a repo
    root reachable from it -- would come back with paths that do not exist. The
    tool's answer must not depend on where the shell happens to be.
    """
    import subprocess

    code = (
        "from pharmacogenomic_record import cli\n"
        "assert cli.PAIRS_PATH.is_file(), cli.PAIRS_PATH\n"
        "assert cli.POSITIONS_PATH.is_file(), cli.POSITIONS_PATH\n"
        "print('ok')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"


# --------------------------------------------------------------------------
# The ingest half, which runs without Docker.
# --------------------------------------------------------------------------


def test_ingest_to_calls_reports_uncovered_genes_without_docker(tmp_path):
    """The ingest half runs without Docker and must surface coverage."""
    vcf_path, report = ingest_to_calls(
        FIXTURES / "23andme_valid.txt", POSITIONS, tmp_path
    )

    assert vcf_path.is_file()
    # The fixture covers two DPYD and two CYP2C19 positions, so CYP2D6 -- which
    # the fixture says nothing about -- must be fully uncovered.
    assert "CYP2D6" in report.genes_fully_uncovered


def test_ingest_of_a_sparse_array_is_mostly_uncovered(tmp_path):
    """The measured, unflattering numbers the README quotes.

    4 of 1226 positions covered is what a real consumer array yields against
    the 3.4.0 reference. It is pinned here so the README's honesty claim cannot
    drift away from the code.
    """
    _vcf_path, report = ingest_to_calls(
        FIXTURES / "23andme_valid.txt", POSITIONS, tmp_path
    )

    assert len(report.covered_rsids) == 4
    assert len(report.uncovered_rsids) == 1014
    assert report.unjoinable_positions == 208
    # Under the strict rule nothing is fully covered by four positions.
    assert report.genes_fully_covered == frozenset()
    assert report.genes_partially_covered == frozenset({"CYP2C19", "DPYD"})


# --------------------------------------------------------------------------
# The strict coverage rule: every joinable position, or the gene is not called.
# --------------------------------------------------------------------------


def test_gene_with_every_joinable_position_covered_is_eligible_to_be_called(
    tmp_path,
):
    """All of CYP2C19's joinable positions covered -> fully covered, not partial.

    This is the rule that makes `called` reachable at all. If
    genes_partially_covered meant "at least one position covered", CYP2C19
    would land in it here and be forced to indeterminate, and no gene could
    ever be called end to end.
    """
    report = build_vcf(CYP2C19_CALLS, REF, tmp_path / "out.vcf")

    assert "CYP2C19" in report.genes_fully_covered
    assert "CYP2C19" not in report.genes_partially_covered
    assert "CYP2C19" not in report.genes_fully_uncovered


def test_gene_one_position_short_is_partial_not_called(tmp_path):
    """One of two joinable positions covered -> partial, never fully covered."""
    report = build_vcf(CYP2C19_CALLS[:1], REF, tmp_path / "out.vcf")

    assert "CYP2C19" in report.genes_partially_covered
    assert "CYP2C19" not in report.genes_fully_covered
    assert "CYP2C19" not in report.genes_fully_uncovered


def test_a_gene_with_no_joinable_position_is_uncovered_not_vacuously_covered(
    tmp_path,
):
    """"Every joinable position covered" is vacuously true of an empty set.

    CFTR in REF is reachable only through an rsID-less position, so its joinable
    set is empty and `set() <= covered` holds. Ranked wrongly -- fully_covered
    tested against every gene rather than only genes with some coverage -- a
    gene the array said nothing whatsoever about would come back eligible to be
    called, which is the worst possible direction for this rule to fail.
    """
    report = build_vcf(CYP2C19_CALLS, REF, tmp_path / "out.vcf")

    assert "CFTR" in report.genes_fully_uncovered
    assert "CFTR" not in report.genes_fully_covered
    assert "CFTR" not in report.genes_partially_covered


def test_the_three_gene_sets_partition_every_gene(tmp_path):
    """Disjoint and exhaustive: no gene may fall through or be double-counted."""
    report = build_vcf(CYP2C19_CALLS[:1], REF, tmp_path / "out.vcf")

    all_genes = {p.gene for p in REF if p.gene is not None}
    sets = (
        report.genes_fully_covered,
        report.genes_partially_covered,
        report.genes_fully_uncovered,
    )

    assert set().union(*sets) == all_genes
    assert sum(len(s) for s in sets) == len(all_genes)


def test_a_fully_covered_gene_reaches_called_end_to_end(tmp_path):
    """Coverage wiring: full coverage lets PharmCAT's call through."""
    report = build_vcf(CYP2C19_CALLS, REF, tmp_path / "out.vcf")

    calls = {c.gene: c for c in calls_from_phenotype(PHENOTYPE_SAMPLE, report)}

    assert calls["CYP2C19"].coverage == CALLED
    assert calls["CYP2C19"].diplotype == "*1/*2"


def test_a_gene_one_position_short_is_recorded_indeterminate(tmp_path):
    """The same PharmCAT output, one position short, must not become a call.

    PharmCAT assumes reference at unobserved positions, so it reports the same
    confident "*1/*2" either way. Coverage is the only thing that knows better.
    """
    report = build_vcf(CYP2C19_CALLS[:1], REF, tmp_path / "out.vcf")

    calls = {c.gene: c for c in calls_from_phenotype(PHENOTYPE_SAMPLE, report)}

    assert calls["CYP2C19"].coverage == INDETERMINATE
    assert calls["CYP2C19"].diplotype is None


def test_a_gene_with_no_covered_position_is_recorded_not_covered(tmp_path):
    report = build_vcf(CYP2C19_CALLS, REF, tmp_path / "out.vcf")

    calls = {c.gene: c for c in calls_from_phenotype(PHENOTYPE_SAMPLE, report)}

    assert calls["DPYD"].coverage == NOT_COVERED


def test_cmd_ingest_prints_coverage_and_surfaces_a_missing_docker(
    store, tmp_path, capsys
):
    """Coverage is reported before PharmCAT is attempted, and the failure raises.

    Docker is absent in this environment, so this exercises the ingest half and
    proves the PharmCAT failure is surfaced rather than turned into an empty
    record.
    """
    with pytest.raises(PharmcatError):
        cmd_ingest(
            store,
            FIXTURES / "23andme_valid.txt",
            "s1",
            tmp_path,
            positions_path=POSITIONS,
        )

    out = capsys.readouterr().out
    assert "covered positions: 4" in out
    assert "uncovered positions: 1014" in out
    assert "genes fully covered (eligible to be called): 0" in out
    assert store.history("s1") == []


def test_error_line_is_ordered_after_the_output_it_follows(tmp_path):
    """Redirected stdout is block-buffered; stderr is not.

    Without an explicit flush the error line lands ABOVE the coverage summary in
    a redirected log, making the failure look like it happened before the work
    that actually preceded it. Run as a subprocess because that is the only way
    to get real buffering behaviour -- capsys replaces the streams.
    """
    import subprocess

    log = tmp_path / "all.log"
    with log.open("w") as handle:
        code = subprocess.call(
            [
                sys.executable, "-m", "pharmacogenomic_record.cli",
                "--db", str(tmp_path / "records.db"),
                "ingest", str(FIXTURES / "23andme_valid.txt"),
                "--subject", "s1",
                "--workdir", str(tmp_path / "work"),
                "--positions", str(POSITIONS),
            ],
            stdout=handle,
            stderr=subprocess.STDOUT,
        )

    lines = log.read_text().splitlines()

    assert code == EXIT_ERROR
    assert lines[0].startswith("wrote ")
    assert lines[-1].startswith("error: ")


def test_main_ingest_exits_nonzero_when_pharmcat_cannot_run(tmp_path, capsys):
    exit_code = main(
        [
            "--db", str(tmp_path / "records.db"),
            "ingest", str(FIXTURES / "23andme_valid.txt"),
            "--subject", "s1",
            "--workdir", str(tmp_path),
            "--positions", str(POSITIONS),
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == EXIT_ERROR
    assert "docker" in captured.err.lower()
    assert_no_reassurance(captured.err)


# --------------------------------------------------------------------------
# The three rendered outcomes.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("subject_id", ["s1", HOSTILE_SUBJECT_ID])
def test_query_prints_cannot_assess_prominently(store, subject_id, capsys):
    store.append(
        subject_id,
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )

    outcome = cmd_query(store, subject_id, "codeine", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert outcome == CANNOT_ASSESS
    assert "CANNOT ASSESS" in out
    assert "no interaction" not in out.lower()
    assert_no_reassurance(out)
    # A cannot_assess answer still cites where the guidance is: the gene being
    # unknown is not a reason to hide that CPIC publishes for this pair.
    assert "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/" in out
    # The subject id is caller text and must not be echoed into the answer.
    assert subject_id not in out


def test_query_prints_guidance_with_citation(store, capsys):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    outcome = cmd_query(store, "s1", "clopidogrel", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert outcome == GUIDANCE_FOUND
    assert "CYP2C19-clopidogrel" in out
    assert "https://" in out
    # A citation, never a recommendation.
    for directive in ("mg", "you should", "dose of", "recommend"):
        assert directive not in out.lower()


def test_query_prints_guidance_for_a_called_gene_with_no_phenotype(store, capsys):
    """VKORC1 has no metabolizer phenotype; coverage alone decides assessability."""
    store.append(
        "s1",
        [
            GeneCall("VKORC1", "rs9923231 variant", None, CALLED),
            GeneCall("CYP2C9", "*1/*1", "Normal Metabolizer", CALLED),
        ],
        guideline_version="cpic-2026-07",
    )

    outcome = cmd_query(store, "s1", "warfarin", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert outcome == GUIDANCE_FOUND
    assert "CANNOT ASSESS" not in out
    assert "VKORC1-warfarin" in out


def test_query_prints_an_explicit_negative_for_a_drug_with_no_cpic_pair(
    store, capsys
):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    cmd_query(store, "s1", "ibuprofen", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert "NO CPIC PAIR" in out
    assert "ibuprofen" in out
    # Silence is the failure mode with consequences.
    assert out.strip()


def test_query_for_an_unknown_subject_cannot_assess(store, capsys):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    outcome = cmd_query(store, "nobody", "clopidogrel", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert outcome == CANNOT_ASSESS
    assert "CANNOT ASSESS" in out
    assert_no_reassurance(out)


def test_query_overall_verdict_is_the_least_reassuring_component(store, capsys):
    """Warfarin has two pairs. One unknown gene invalidates the whole answer."""
    store.append(
        "s1",
        [
            GeneCall("CYP2C9", "*1/*1", "Normal Metabolizer", CALLED),
            GeneCall("VKORC1", None, None, NOT_COVERED),
        ],
        guideline_version="cpic-2026-07",
    )

    outcome = cmd_query(store, "s1", "warfarin", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert outcome == CANNOT_ASSESS
    assert "OVERALL: CANNOT ASSESS" in out
    # Both components are still reported; the summary does not replace them.
    assert "CYP2C9-warfarin" in out
    assert "VKORC1-warfarin" in out


def test_query_output_says_it_is_not_a_medical_device(store, capsys):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    cmd_query(store, "s1", "clopidogrel", pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert "not a medical device" in out.lower()


# --------------------------------------------------------------------------
# Exit codes.
# --------------------------------------------------------------------------


def _seed(tmp_path, calls):
    db = tmp_path / "records.db"
    RecordStore(db).append("s1", calls, guideline_version="cpic-2026-07")
    return db


def test_exit_code_distinguishes_cannot_assess_from_an_answer(tmp_path, capsys):
    db = _seed(tmp_path, [GeneCall("CYP2D6", None, None, NOT_COVERED)])

    code = main(["--db", str(db), "query", "codeine", "--subject", "s1",
                 "--pairs", str(PAIRS)])
    capsys.readouterr()

    assert code == EXIT_CANNOT_ASSESS
    assert code != EXIT_OK


def test_exit_code_is_zero_for_an_assessed_query(tmp_path, capsys):
    db = _seed(
        tmp_path,
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
    )

    assert (
        main(["--db", str(db), "query", "clopidogrel", "--subject", "s1",
              "--pairs", str(PAIRS)])
        == EXIT_OK
    )
    capsys.readouterr()


def test_exit_code_is_zero_for_a_drug_absent_from_the_table(tmp_path, capsys):
    """An explicit statement about the table is an answer, not an unknown."""
    db = _seed(
        tmp_path,
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
    )

    assert (
        main(["--db", str(db), "query", "ibuprofen", "--subject", "s1",
              "--pairs", str(PAIRS)])
        == EXIT_OK
    )
    capsys.readouterr()


# --------------------------------------------------------------------------
# Failures are surfaced, never rendered as "nothing to report".
# --------------------------------------------------------------------------


def test_an_unloadable_pair_table_is_surfaced_not_reported_as_no_guidance(
    tmp_path, capsys
):
    bad = tmp_path / "pairs.json"
    bad.write_text(json.dumps([{"gene": "CYP2C19"}]), encoding="utf-8")
    db = _seed(
        tmp_path,
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
    )

    code = main(["--db", str(db), "query", "clopidogrel", "--subject", "s1",
                 "--pairs", str(bad)])
    captured = capsys.readouterr()

    assert code == EXIT_ERROR
    assert "missing required field" in captured.err
    assert "NO CPIC PAIR" not in captured.out
    assert_no_reassurance(captured.out)
    assert_no_reassurance(captured.err)


def test_a_missing_pair_table_is_surfaced(tmp_path, capsys):
    db = _seed(
        tmp_path,
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
    )

    code = main(["--db", str(db), "query", "clopidogrel", "--subject", "s1",
                 "--pairs", str(tmp_path / "absent.json")])
    captured = capsys.readouterr()

    assert code == EXIT_ERROR
    assert "could not read the gene-drug pair table" in captured.err


def test_a_lowercase_gene_is_rejected_at_write_time_and_surfaced(
    tmp_path, capsys, monkeypatch
):
    """PharmCAT 3.4.0 emits uppercase; if that ever changes, it must not pass.

    The store's CHECK is what rejects it. What is tested here is that the CLI
    lets that rejection out as a non-zero exit with the message on stderr,
    rather than swallowing it and reporting a successful ingest.
    """
    import pharmacogenomic_record.cli as cli

    monkeypatch.setattr(cli, "run_pharmcat", lambda vcf, workdir: PHENOTYPE_SAMPLE)
    monkeypatch.setattr(
        cli,
        "calls_from_phenotype",
        lambda phenotype_json, report: [
            GeneCall("cyp2c19", "*1/*2", "Intermediate Metabolizer", CALLED)
        ],
    )

    code = main(
        [
            "--db", str(tmp_path / "records.db"),
            "ingest", str(FIXTURES / "23andme_valid.txt"),
            "--subject", "s1",
            "--workdir", str(tmp_path / "work"),
            "--positions", str(POSITIONS),
        ]
    )
    captured = capsys.readouterr()

    assert code == EXIT_ERROR
    assert "CHECK constraint failed" in captured.err
    assert "stored record" not in captured.out
    assert RecordStore(tmp_path / "records.db").history("s1") == []


# --------------------------------------------------------------------------
# Drift, wired with a real collection of pair ids.
# --------------------------------------------------------------------------


def test_drift_reports_affected_subjects_including_uncovered_genes(store, capsys):
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )

    cmd_drift(store, ["CYP2D6-codeine"], pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert "s1" in out
    assert "CYP2D6" in out
    assert "CYP2D6-codeine" in out


def test_drift_marks_each_affected_record_so_a_dropped_list_fails_here_too(
    store, capsys
):
    """A second CLI test on the affected branch.

    A regression that dropped the affected records (`affected = []`) was caught
    by only one CLI test, because the "no stored record" negative it fell back
    to still echoes the pair id every other assertion looked for. This anchors
    on the `[AFFECTED]` marker, which the negative branch never prints, so the
    dropped-records path is caught by more than one test.
    """
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )

    cmd_drift(store, ["CYP2D6-codeine"], pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert "[AFFECTED]" in out
    assert "no stored record" not in out.lower()


def test_drift_prints_an_explicit_negative_when_nobody_is_affected(store, capsys):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    cmd_drift(store, ["CYP2D6-codeine"], pairs_path=PAIRS)
    out = capsys.readouterr().out

    assert "no stored record" in out.lower()
    assert out.strip()


def test_drift_unknown_pair_id_is_not_reported_as_a_genuine_negative(store, capsys):
    """A typo'd pair id must warn, not borrow the shape of "touches nobody".

    `--changed-pair CYP2D6-codiene` (a misspelling absent from the table) and a
    real pair that simply matches no stored record are two different facts: the
    first means "your id matched nothing", the second means "the revision
    touches nobody stored". Printing the same reassuring negative for a typo is
    the same class of danger as `cannot_assess` reading like an all-clear -- the
    user thinks they got an answer about the pair they meant.
    """
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", CALLED)],
        guideline_version="cpic-2026-07",
    )

    # A genuine negative: a real pair, but no stored record holds its gene.
    cmd_drift(store, ["CYP2D6-codeine"], pairs_path=PAIRS)
    genuine = capsys.readouterr().out

    # A typo: the id matches no gene-drug pair in the table at all.
    cmd_drift(store, ["CYP2D6-codiene"], pairs_path=PAIRS)
    typo = capsys.readouterr().out

    # The two negatives must not read identically.
    assert genuine != typo
    # The genuine negative is the "revision touches nobody stored" statement.
    assert "no stored record holds" in genuine.lower()
    # The typo must NOT reuse that statement; that would read as an all-clear
    # about a pair that was never actually checked.
    assert "no stored record holds" not in typo.lower()
    # It names the bad id, warns, and does not read as reassurance.
    assert "CYP2D6-codiene" in typo
    assert_no_reassurance(typo)


def test_drift_warns_on_a_typo_even_alongside_a_real_match(store, capsys):
    """A good id and a bad id together: report the match, still flag the typo.

    The warning about the unknown id must not be swallowed just because another
    id in the same invocation did match something.
    """
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, NOT_COVERED)],
        guideline_version="cpic-2026-07",
    )

    cmd_drift(store, ["CYP2D6-codeine", "CYP2D6-codiene"], pairs_path=PAIRS)
    out = capsys.readouterr().out

    # The real pair still produces its affected line...
    assert "[AFFECTED]" in out
    assert "s1" in out
    # ...and the typo is still flagged, naming only the unknown id.
    assert "CYP2D6-codiene" in out
    assert_no_reassurance(out)


def test_main_drift_passes_a_collection_not_a_bare_string(tmp_path, capsys):
    """`--changed-pair` must reach drift as a collection; a str raises TypeError."""
    db = _seed(tmp_path, [GeneCall("CYP2D6", None, None, NOT_COVERED)])

    code = main(["--db", str(db), "drift", "--changed-pair", "CYP2D6-codeine",
                 "--pairs", str(PAIRS)])
    out = capsys.readouterr().out

    assert code == EXIT_OK
    assert "CYP2D6-codeine" in out
