import json
from pathlib import Path

import pytest

from pharmacogenomic_record.caller import GeneCall, PharmcatError, parse_phenotype_json, run_pharmcat

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


def test_indeterminate_phenotype_is_not_called(tmp_path):
    """PharmCAT saying 'Indeterminate' is not a normal-metabolizer result.

    Both alleles are named here, so the structural no-call branch cannot fire:
    this reaches the phenotype-marker branch and nothing else. The indeterminate
    state must also null both payload fields, otherwise a leaked label reads as
    a call downstream.
    """
    path = write_payload(
        tmp_path,
        {"TPMT": {"gene": "TPMT", "diplotypes": [
            diplotype("*1/*3A", "*1", "*3A", ["Indeterminate"])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


def test_committed_fixture_null_allele_gene_is_indeterminate():
    """The fixture's TPMT has null alleles, so the STRUCTURAL branch fires.

    Kept separate from the phenotype-marker test above: this one pins the
    committed realistic sample, where the no-call is expressed structurally.
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


def test_same_label_conflicting_phenotypes_is_indeterminate(tmp_path):
    """Agreeing on the label is not agreeing on the call.

    Two candidates labelled *1/*2 carrying Poor and Normal Metabolizer are
    clinically opposite. Keying ambiguity on the label alone kept whichever
    happened to be first and silently discarded the other.
    """
    path = write_payload(
        tmp_path,
        {
            "CYP2C19": {
                "gene": "CYP2C19",
                "diplotypes": [
                    diplotype("*1/*2", "*1", "*2", ["Poor Metabolizer"]),
                    diplotype("*1/*2", "*1", "*2", ["Normal Metabolizer"]),
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


def test_same_label_conflicting_allele_names_is_indeterminate(tmp_path):
    """Two candidates can share a label while naming different alleles."""
    path = write_payload(
        tmp_path,
        {
            "CYP2C19": {
                "gene": "CYP2C19",
                "diplotypes": [
                    diplotype("*1/*2", "*1", "*2", ["Intermediate Metabolizer"]),
                    diplotype("*1/*2", "*2", "*1", ["Intermediate Metabolizer"]),
                ],
            }
        },
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None


def test_a_no_call_candidate_anywhere_blocks_the_call(tmp_path):
    """The structural no-call check applies to EVERY candidate, not just the first.

    A named candidate followed by an all-null candidate under the same label
    used to be reported as a confident call, discarding the no-call entirely.
    """
    path = write_payload(
        tmp_path,
        {
            "CYP2C19": {
                "gene": "CYP2C19",
                "diplotypes": [
                    diplotype("*1/*2", "*1", "*2", ["Normal Metabolizer"]),
                    {
                        "allele1": None,
                        "allele2": None,
                        "label": "*1/*2",
                        "phenotypes": ["Indeterminate"],
                    },
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


def test_non_string_label_on_a_later_candidate_raises(tmp_path):
    """Unparseable output must raise wherever it appears, not degrade quietly.

    Only the first candidate's label used to be type-checked, so a non-string
    label further down came out as a plain 'indeterminate'.
    """
    path = write_payload(
        tmp_path,
        {
            "CYP2C19": {
                "gene": "CYP2C19",
                "diplotypes": [
                    diplotype("*1/*2", "*1", "*2", ["Normal Metabolizer"]),
                    {"label": 7},
                ],
            }
        },
    )
    with pytest.raises(PharmcatError, match="non-string label"):
        parse_phenotype_json(
            path, uncovered_genes=set(), partially_covered_genes=set()
        )


def test_named_alleles_with_a_null_label_is_indeterminate(tmp_path):
    """A call needs something to report. Named alleles alone are not enough."""
    path = write_payload(
        tmp_path,
        {"CYP2C19": {"gene": "CYP2C19", "diplotypes": [
            diplotype(None, "*1", "*2", ["Normal Metabolizer"])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "indeterminate"
    assert call.diplotype is None
    assert call.phenotype is None


@pytest.mark.parametrize(
    "phenotype",
    ["No Call", "Not Available", "No Data", "NA", "Undetermined", "Not Assigned"],
)
def test_additional_indeterminate_markers_downgrade_a_complete_diplotype(
    tmp_path, phenotype
):
    """These real PharmCAT strings mean 'no result' and must not read as calls."""
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


@pytest.mark.parametrize(
    "phenotype",
    [
        "Normal Metabolizer",
        "Intermediate Metabolizer",
        "Poor Metabolizer",
        "Ultrarapid Metabolizer",
        "Normal Function",
        "No Function",
        "Decreased Function",
        "Increased Function",
        "Possible Decreased Function",
        "Deficient",
        "Deficient with CNSHA",
        "Variable",
        "Uncertain Susceptibility",
        "Malignant Hyperthermia Susceptibility",
    ],
)
def test_legitimate_cpic_phenotypes_are_never_downgraded(tmp_path, phenotype):
    """The marker list must not swallow real results.

    'No Function' vs the 'No Data'/'No Call' markers is the sharp edge: they
    share the leading 'no' token, and only whole-token run matching keeps a
    no-function allele -- clinically actionable -- out of the marker list.
    """
    path = write_payload(
        tmp_path,
        {"CYP2C9": {"gene": "CYP2C9", "diplotypes": [
            diplotype("*2/*3", "*2", "*3", [phenotype])
        ]}},
    )
    (call,) = parse_phenotype_json(
        path, uncovered_genes=set(), partially_covered_genes=set()
    )

    assert call.coverage == "called"
    assert call.diplotype == "*2/*3"
    assert call.phenotype == phenotype


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


def test_the_three_coverage_states_do_not_collapse():
    """Assert the exact state per gene, not a subset of the allowed values.

    A subset assertion passes even if every gene collapsed to one state, which
    is precisely the failure this software must never have.
    """
    calls = parse_phenotype_json(
        SAMPLE, uncovered_genes={"CYP2D6"}, partially_covered_genes={"NAT2"}
    )
    states = {c.gene: c.coverage for c in calls}

    assert states == {
        "CYP2C19": "called",       # fully covered, unambiguous
        "DPYD": "called",          # fully covered, unambiguous
        "TPMT": "indeterminate",   # null alleles in the fixture
        "CYP2D6": "not_covered",   # declared uncovered, absent from output
        "NAT2": "indeterminate",   # declared partially covered
    }
    # All three states are genuinely present, so none has swallowed the others.
    assert set(states.values()) == {"called", "not_covered", "indeterminate"}


def test_malformed_json_raises(tmp_path):
    bad = tmp_path / "pharmcat_malformed.json"
    bad.write_text("{not json")
    with pytest.raises(PharmcatError, match="parse"):
        parse_phenotype_json(
            bad, uncovered_genes=set(), partially_covered_genes=set()
        )


@pytest.mark.parametrize(
    "payload",
    [
        {"totally": "wrong"},   # key absent -> phenotypes is None
        {"phenotypes": []},     # present but a list, a distinct wrong shape
    ],
    ids=["key-absent", "list-not-object"],
)
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
    from pharmacogenomic_record.caller import _vcf_basename

    assert _vcf_basename(Path("/tmp") / name) == expected
