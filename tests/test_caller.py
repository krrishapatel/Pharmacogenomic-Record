import json
from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall, PharmcatError, parse_phenotype_json, run_pharmcat

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SAMPLE = FIXTURES / "pharmcat_phenotype_sample.json"


def write_payload(tmp_path, phenotypes, name="phenotype.json"):
    """Write an inline PharmCAT-shaped payload so the committed fixture stays real."""
    path = tmp_path / name
    path.write_text(json.dumps({"phenotypes": phenotypes}))
    return path


def diplotype(label, allele1, allele2, phenotypes):
    return {
        "allele1": {"name": allele1},
        "allele2": {"name": allele2},
        "label": label,
        "phenotypes": phenotypes,
    }


def test_parses_called_genes():
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes=set(), partially_covered_genes=set()
    )
    by_gene = {c.gene: c for c in calls}

    assert by_gene["CYP2C19"] == GeneCall(
        gene="CYP2C19",
        diplotype="*1/*2",
        phenotype="Intermediate Metabolizer",
        coverage="called",
    )


def test_indeterminate_phenotype_is_not_called():
    """PharmCAT saying 'Indeterminate' is not a normal-metabolizer result.

    The indeterminate state must also null both payload fields, otherwise a
    leaked label reads as a call downstream.
    """
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes=set(), partially_covered_genes=set()
    )
    tpmt = next(c for c in calls if c.gene == "TPMT")

    assert tpmt.coverage == "indeterminate"
    assert tpmt.diplotype is None
    assert tpmt.phenotype is None


def test_null_alleles_are_the_no_call_signal(tmp_path):
    """A null allele is PharmCAT's structural no-call, whatever the label says."""
    path = write_payload(
        tmp_path,
        {
            "TPMT": {
                "gene": "TPMT",
                "diplotypes": [
                    {
                        "allele1": None,
                        "allele2": None,
                        "label": "*1/*1",
                        "phenotypes": ["Normal Metabolizer"],
                    }
                ],
            }
        },
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


@pytest.mark.parametrize(
    "allele2",
    [
        {"name": "  "},  # blank name
        {"name": None},  # explicitly null name
        {},  # name key missing entirely
    ],
    ids=["blank", "null", "missing"],
)
def test_unnamed_allele_is_a_no_call(tmp_path, allele2):
    """A null/missing/blank allele name is a no-call, not a parse error."""
    path = write_payload(
        tmp_path,
        {
            "TPMT": {
                "gene": "TPMT",
                "diplotypes": [
                    {
                        "allele1": {"name": "*1"},
                        "allele2": allele2,
                        "label": "*1/*1",
                        "phenotypes": ["Normal Metabolizer"],
                    }
                ],
            }
        },
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )
    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


@pytest.mark.parametrize(
    "phenotype",
    ["Indeterminate", "Indeterminate Metabolizer", "N/A", "Unknown", "No Result"],
)
def test_indeterminate_phenotype_string_downgrades_a_complete_diplotype(
    tmp_path, phenotype
):
    """Marker matching still applies -- to the phenotype, on whole tokens.

    The alleles here are fully named, so only the phenotype string can trigger
    this. 'N/A' must match as the two tokens n + a, not as a raw substring.
    """
    path = write_payload(
        tmp_path,
        {"CYP2C19": {"gene": "CYP2C19", "diplotypes": [
            diplotype("*1/*2", "*1", "*2", [phenotype])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


def test_multiple_candidate_diplotypes_are_indeterminate(tmp_path):
    """Unphased data consistent with several allele combinations is ambiguity.

    Reporting only diplotypes[0] would present *1/*4 Intermediate as THE call
    when *2/*41 Normal is equally consistent.
    """
    path = write_payload(
        tmp_path,
        {
            "CYP2D6": {
                "gene": "CYP2D6",
                "diplotypes": [
                    diplotype("*1/*4", "*1", "*4", ["Intermediate Metabolizer"]),
                    diplotype("*2/*41", "*2", "*41", ["Normal Metabolizer"]),
                ],
            }
        },
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


def test_identical_candidate_labels_still_called(tmp_path):
    """Duplicate labels are not ambiguity -- do not over-trigger on list length."""
    path = write_payload(
        tmp_path,
        {
            "CYP2C19": {
                "gene": "CYP2C19",
                "diplotypes": [
                    diplotype("*1/*2", "*1", "*2", ["Intermediate Metabolizer"]),
                    diplotype("*1/*2", "*1", "*2", ["Intermediate Metabolizer"]),
                ],
            }
        },
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "called"
    assert call.diplotype == "*1/*2"
    assert call.phenotype == "Intermediate Metabolizer"


@pytest.mark.parametrize(
    ("allele1", "allele2", "label"),
    [
        ("Canton", "Aures", "Canton/Aures"),
        ("Mediterranean", "Asahi", "Mediterranean/Asahi"),
        ("Viangchan", "A-202A_376G", "Viangchan/A-202A_376G"),
    ],
)
def test_g6pd_place_named_alleles_are_called(tmp_path, allele1, allele2, label):
    """Regression: the 'n/a' marker used to straddle the label's separator.

    'n/a' in 'canton/aures' is True, so real G6PD-deficient results were being
    silently discarded as indeterminate. The label must never be pattern-matched.
    """
    path = write_payload(
        tmp_path,
        {"G6PD": {"gene": "G6PD", "diplotypes": [
            diplotype(label, allele1, allele2, ["Deficient"])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "called"
    assert call.diplotype == label
    assert call.phenotype == "Deficient"


def test_gene_without_a_phenotype_is_still_called(tmp_path):
    """F2, F5, VKORC1, CFTR, IFNL3 and ABCG2 have no metabolizer phenotype.

    Absence of a metabolizer phenotype is not absence of a call.
    """
    path = write_payload(
        tmp_path,
        {"F5": {"gene": "F5", "diplotypes": [
            {
                "allele1": {"name": "rs6025 reference (C)"},
                "allele2": {"name": "rs6025 reference (C)"},
                "label": "rs6025 reference (C)/rs6025 reference (C)",
            }
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "called"
    assert call.diplotype == "rs6025 reference (C)/rs6025 reference (C)"
    assert call.phenotype is None


def test_no_function_phenotype_is_not_indeterminate(tmp_path):
    """'No Function' must not trip the 'no result' marker."""
    path = write_payload(
        tmp_path,
        {"CYP2C9": {"gene": "CYP2C9", "diplotypes": [
            diplotype("*3/*3", "*3", "*3", ["No Function"])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )
    assert call.coverage == "called"


def test_empty_diplotype_list_is_indeterminate(tmp_path):
    path = write_payload(tmp_path, {"NAT2": {"gene": "NAT2", "diplotypes": []}})
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )
    assert call.coverage == "indeterminate"


def test_uncovered_genes_are_marked_not_covered():
    """A gene absent from the array must never be reported as called."""
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes={"CYP2D6", "DPYD"}, partially_covered_genes=set()
    )
    by_gene = {c.gene: c for c in calls}

    assert by_gene["CYP2D6"].coverage == "not_covered"
    assert by_gene["CYP2D6"].phenotype is None
    assert by_gene["CYP2D6"].diplotype is None
    # DPYD appears in the JSON but the array did not cover it -- coverage wins.
    assert by_gene["DPYD"].coverage == "not_covered"
    assert by_gene["DPYD"].phenotype is None


def test_partially_covered_gene_is_indeterminate():
    """PharmCAT assumes reference at unobserved positions, so its confident
    Normal Metabolizer for a barely-measured gene must not be trusted."""
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes=set(), partially_covered_genes={"DPYD"}
    )
    dpyd = next(c for c in calls if c.gene == "DPYD")

    assert dpyd.coverage == "indeterminate"
    assert dpyd.diplotype is None
    assert dpyd.phenotype is None


def test_partially_covered_gene_absent_from_pharmcat_output():
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes=set(), partially_covered_genes={"NUDT15"}
    )
    nudt15 = next(c for c in calls if c.gene == "NUDT15")
    assert nudt15 == GeneCall("NUDT15", None, None, "indeterminate")


def test_uncovered_beats_partially_covered():
    """Precedence is strict: uncovered > partial > PharmCAT output."""
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes={"DPYD"}, partially_covered_genes={"DPYD"}
    )
    dpyd = [c for c in calls if c.gene == "DPYD"]

    assert len(dpyd) == 1
    assert dpyd[0].coverage == "not_covered"


def test_uncovered_beats_partially_covered_when_absent_from_output():
    """Precedence must hold for the synthesized rows too, and emit exactly one."""
    calls = parse_phenotype_json(
        SAMPLE,
        uncovered_genes={"NUDT15"},
        partially_covered_genes={"NUDT15"},
    )
    nudt15 = [c for c in calls if c.gene == "NUDT15"]

    assert nudt15 == [GeneCall("NUDT15", None, None, "not_covered")]


def test_no_gene_is_reported_twice():
    calls = parse_phenotype_json(
        SAMPLE,
        uncovered_genes={"CYP2D6", "DPYD"},
        partially_covered_genes={"CYP2D6", "DPYD", "NAT2", "CYP2C19"},
    )
    genes = [c.gene for c in calls]
    assert len(genes) == len(set(genes))


def test_every_coverage_value_is_one_of_three_states():
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes={"CYP2D6"}, partially_covered_genes={"NAT2"}
    )
    assert {c.coverage for c in calls} <= {"called", "not_covered", "indeterminate"}


def test_malformed_json_raises(tmp_path):
    bad = tmp_path / "pharmcat_malformed.json"
    bad.write_text("{not json")
    with pytest.raises(PharmcatError, match="parse"):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


@pytest.mark.parametrize("payload", [{"totally": "wrong"}, {}, {"phenotypes": None}])
def test_wrong_shaped_payload_raises(tmp_path, payload):
    """Returning [] here is indistinguishable from 'PharmCAT called nothing'."""
    bad = tmp_path / "wrong.json"
    bad.write_text(json.dumps(payload))
    with pytest.raises(PharmcatError, match="phenotypes"):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


def test_non_dict_json_raises(tmp_path):
    bad = tmp_path / "list.json"
    bad.write_text("[1, 2, 3]")
    with pytest.raises(PharmcatError, match="not a JSON object"):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


def test_non_utf8_output_raises(tmp_path):
    bad = tmp_path / "binary.json"
    bad.write_bytes(b"\xff\xfe\x00not utf-8")
    with pytest.raises(PharmcatError, match="parse"):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


@pytest.mark.parametrize(
    "entry",
    [
        "oops",
        {"diplotypes": "oops"},
        {"diplotypes": ["oops"]},
        {"diplotypes": [{"allele1": "oops", "allele2": {"name": "*1"}}]},
        {"diplotypes": [{"allele1": {"name": 7}, "allele2": {"name": "*1"}}]},
    ],
)
def test_unparseable_gene_entry_raises_pharmcat_error(tmp_path, entry):
    """Every unparseable-output path must surface as PharmcatError."""
    bad = tmp_path / "bad_entry.json"
    bad.write_text(json.dumps({"phenotypes": {"TPMT": entry}}))
    with pytest.raises(PharmcatError):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


def test_empty_result_with_no_coverage_information_raises(tmp_path):
    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"phenotypes": {}}))
    with pytest.raises(PharmcatError, match="nothing to record"):
        parse_phenotype_json(
            empty, uncovered_genes=set(), partially_covered_genes=set()
        )


def test_empty_pharmcat_output_with_coverage_sets_does_not_raise(tmp_path):
    """All-not_covered / all-indeterminate rows are a legitimate result."""
    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"phenotypes": {}}))
    calls = parse_phenotype_json(
        empty, uncovered_genes={"CYP2D6"}, partially_covered_genes={"CYP2C9"}
    )

    assert calls == [
        GeneCall("CYP2D6", None, None, "not_covered"),
        GeneCall("CYP2C9", None, None, "indeterminate"),
    ]


def test_coverage_arguments_accept_any_iterable(tmp_path):
    calls = parse_phenotype_json(
        SAMPLE,
        uncovered_genes=(g for g in ["CYP2D6"]),
        partially_covered_genes=["NAT2"],
    )
    by_gene = {c.gene: c for c in calls}
    assert by_gene["CYP2D6"].coverage == "not_covered"
    assert by_gene["NAT2"].coverage == "indeterminate"


def test_coverage_arguments_are_required():
    """A coverage guard that can be omitted is a guard that will be omitted."""
    with pytest.raises(TypeError):
        parse_phenotype_json(SAMPLE)  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        parse_phenotype_json(SAMPLE, set())  # type: ignore[call-arg]


def test_run_pharmcat_raises_when_docker_missing(tmp_path, monkeypatch):
    """No record may be written when the caller cannot run."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    monkeypatch.setenv("PATH", str(tmp_path))  # hide docker

    with pytest.raises(PharmcatError):
        run_pharmcat(vcf, tmp_path)


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("sample.vcf", "sample"),
        ("sample.vcf.gz", "sample"),
        ("sample.with.dots.vcf.gz", "sample.with.dots"),
        ("sample", "sample"),
    ],
)
def test_vcf_basename_strips_double_extension(name, expected):
    from pgxrecord.caller import _vcf_basename

    assert _vcf_basename(Path("/tmp") / name) == expected
