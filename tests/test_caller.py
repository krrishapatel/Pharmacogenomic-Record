from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall, PharmcatError, parse_phenotype_json, run_pharmcat

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SAMPLE = FIXTURES / "pharmcat_phenotype_sample.json"


def test_parses_called_genes():
    calls = parse_phenotype_json(SAMPLE, uncovered_genes=set())
    by_gene = {c.gene: c for c in calls}

    assert by_gene["CYP2C19"] == GeneCall(
        gene="CYP2C19",
        diplotype="*1/*2",
        phenotype="Intermediate Metabolizer",
        coverage="called",
    )


def test_indeterminate_phenotype_is_not_called():
    """PharmCAT saying 'Indeterminate' is not a normal-metabolizer result."""
    calls = parse_phenotype_json(SAMPLE, uncovered_genes=set())
    tpmt = next(c for c in calls if c.gene == "TPMT")

    assert tpmt.coverage == "indeterminate"
    assert tpmt.phenotype != "Normal Metabolizer"


def test_uncovered_genes_are_marked_not_covered():
    """A gene absent from the array must never be reported as called."""
    calls = parse_phenotype_json(SAMPLE, uncovered_genes={"CYP2D6", "DPYD"})
    by_gene = {c.gene: c for c in calls}

    assert by_gene["CYP2D6"].coverage == "not_covered"
    assert by_gene["CYP2D6"].phenotype is None
    assert by_gene["CYP2D6"].diplotype is None
    # DPYD appears in the JSON but the array did not cover it -- coverage wins.
    assert by_gene["DPYD"].coverage == "not_covered"
    assert by_gene["DPYD"].phenotype is None


def test_every_coverage_value_is_one_of_three_states():
    calls = parse_phenotype_json(SAMPLE, uncovered_genes={"CYP2D6"})
    assert {c.coverage for c in calls} <= {"called", "not_covered", "indeterminate"}


def test_malformed_json_raises(tmp_path):
    bad = tmp_path / "pharmcat_malformed.json"
    bad.write_text("{not json")
    with pytest.raises(PharmcatError, match="parse"):
        parse_phenotype_json(bad, uncovered_genes=set())


def test_run_pharmcat_raises_when_docker_missing(tmp_path, monkeypatch):
    """No record may be written when the caller cannot run."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    monkeypatch.setenv("PATH", str(tmp_path))  # hide docker

    with pytest.raises(PharmcatError):
        run_pharmcat(vcf, tmp_path)
