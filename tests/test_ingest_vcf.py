
import pytest

from pharmacogenomic_record.ingest.raw import RawCall
from pharmacogenomic_record.ingest.vcf import build_vcf
from pharmacogenomic_record.positions import ReferencePosition

REF = [
    ReferencePosition(
        chrom="chr1", pos=100, rsid="rs1", ref="G", alt=("T",), gene="DPYD"
    ),
    ReferencePosition(
        chrom="chr1", pos=200, rsid="rs2", ref="C", alt=("A",), gene="DPYD"
    ),
    ReferencePosition(
        chrom="chr10", pos=300, rsid="rs3", ref="C", alt=("T",), gene="CYP2C19"
    ),
    ReferencePosition(
        chrom="chr22", pos=400, rsid=None, ref="A", alt=("G",), gene="CYP2D6"
    ),
]


def test_writes_vcf_with_matched_positions(tmp_path):
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=999, genotype="GT")]

    build_vcf(calls, REF, out)
    text = out.read_text()

    assert text.startswith("##fileformat=VCFv4.2")
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" in text
    # Every contig we emit must be declared in the header.
    assert '##contig=<ID=chr1,assembly=GRCh38.p14,species="Homo sapiens">' in text
    # Uses the GRCh38 coordinate from the reference, NOT the raw file's 999.
    assert "chr1\t100\trs1\tG\tT\t.\tPASS\t.\tGT\t0/1" in text
    assert "999" not in text


def test_contigs_and_rows_are_in_natural_chromosome_order(tmp_path):
    """chr10 must not sort before chr2. Plain string sort gets this wrong."""
    out = tmp_path / "out.vcf"
    calls = [
        RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG"),
        RawCall(rsid="rs3", chrom="10", pos=300, genotype="CC"),
    ]

    build_vcf(calls, REF, out)
    lines = out.read_text().splitlines()

    contigs = [line for line in lines if line.startswith("##contig")]
    assert "ID=chr1," in contigs[0]
    assert "ID=chr10," in contigs[1]
    assert "ID=chr22," in contigs[2]

    data = [line for line in lines if not line.startswith("#")]
    assert [line.split("\t")[0] for line in data] == ["chr1", "chr10"]


def test_genotype_translation():
    """Raw allele letters become VCF numeric genotypes against the ref allele."""
    from pharmacogenomic_record.ingest.vcf import translate_genotype

    ref = REF[0]  # ref=G alt=T
    assert translate_genotype("GG", ref) == "0/0"
    assert translate_genotype("GT", ref) == "0/1"
    assert translate_genotype("TG", ref) == "0/1"
    assert translate_genotype("TT", ref) == "1/1"
    assert translate_genotype("AA", ref) is None  # allele not in ref/alt


def test_coverage_report_distinguishes_partial_from_absent(tmp_path):
    out = tmp_path / "out.vcf"
    # rs1 present (DPYD partial), rs2 absent, rs3 absent (CYP2C19 fully absent)
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    assert report.covered_rsids == {"rs1"}
    assert report.uncovered_rsids == {"rs2", "rs3"}
    assert "CYP2C19" in report.genes_fully_uncovered
    assert "DPYD" in report.genes_partially_covered
    assert "DPYD" not in report.genes_fully_uncovered


def test_gene_whose_positions_all_lack_rsids_is_fully_uncovered(tmp_path):
    """A gene reachable only via rsID-less positions can never be covered.

    This tests the join mechanism, not a claim about any specific gene. In the
    real 3.4.0 reference NO gene is entirely rsID-less -- CYP2D6, for
    instance, has 146 of 157 positions carrying rsIDs. CYP2D6 is genuinely
    unresolvable from a consumer array, but because it depends on copy-number
    and structural variation, which is a different limitation than the rsID
    join and is not something this function can detect.

    The synthetic CYP2D6 entry in REF models the rsID-less case so the
    mechanism is covered; do not read it as a fact about real CYP2D6 data.
    """
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    assert "CYP2D6" in report.genes_fully_uncovered


def test_unjoinable_rsid_less_positions_are_counted(tmp_path):
    """The 208 rsID-less positions are in neither covered nor uncovered.

    They have no rsID to key on, so they cannot appear in either rsID set --
    but they are real coverage gaps and must be visible, not silently
    dropped. Downstream code must treat "in neither set" as uncovered.
    """
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    # REF has one rsid=None position (the synthetic CYP2D6 entry).
    assert report.unjoinable_positions == 1
    assert report.covered_rsids | report.uncovered_rsids == {"rs1", "rs2", "rs3"}


def test_coverage_report_cannot_be_mutated(tmp_path):
    """frozen=True alone would not stop report.covered_rsids.add(...)."""
    out = tmp_path / "out.vcf"
    report = build_vcf(
        [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")], REF, out
    )

    assert isinstance(report.covered_rsids, frozenset)
    with pytest.raises(AttributeError):
        report.covered_rsids.add("rs99")


def test_partially_covered_gene_is_not_evidence_of_absence(tmp_path):
    """genes_partially_covered means >=1 position, never "fully assessed".

    DPYD has two positions in REF and only one is covered. Treating that as
    a complete assessment would let a missed position read as "no
    interaction found" -- the exact collapse the invariant forbids.
    """
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    assert "DPYD" in report.genes_partially_covered
    assert "rs2" in report.uncovered_rsids  # the other DPYD position


def test_hemizygous_call_is_uncovered_but_reported_as_measured(tmp_path):
    """A single-allele call is not covered, and not silently absent either.

    It writes no VCF row, so it cannot count as covered. But the array did
    report it, and folding that into plain "uncovered" would tell the user their
    array never measured a position it did measure. For a male sample this is
    every joinable G6PD position.
    """
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="X", pos=100, genotype="G")]

    report = build_vcf(calls, REF, out)

    assert "rs1" not in report.covered_rsids
    assert "rs1" in report.uncovered_rsids
    assert "rs1" not in out.read_text()

    assert report.hemizygous_rsids == frozenset({"rs1"})
    assert report.hemizygous_genes == frozenset({"DPYD"})


def test_hemizygous_set_is_empty_for_ordinary_diploid_input(tmp_path):
    """The new sets stay empty rather than shadowing covered positions."""
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    assert report.covered_rsids == frozenset({"rs1"})
    assert report.hemizygous_rsids == frozenset()
    assert report.hemizygous_genes == frozenset()


def test_hemizygous_set_ignores_positions_the_table_does_not_need(tmp_path):
    """A consumer file is full of hemizygous chrY and chrM calls.

    Counting all of them would report a large number that says nothing about
    pharmacogenomic coverage, so only positions the reference table actually
    lists are counted.
    """
    out = tmp_path / "out.vcf"
    calls = [
        RawCall(rsid="rs1", chrom="X", pos=100, genotype="G"),
        RawCall(rsid="rs999999", chrom="Y", pos=1, genotype="A"),
    ]

    report = build_vcf(calls, REF, out)

    assert report.hemizygous_rsids == frozenset({"rs1"})


def test_untranslatable_genotype_counts_as_uncovered(tmp_path):
    """A call whose alleles don't match ref/alt yields no data, not a ref call."""
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="AA")]

    report = build_vcf(calls, REF, out)

    assert "rs1" not in report.covered_rsids
    assert "rs1" in report.uncovered_rsids
    assert "rs1" not in out.read_text()
