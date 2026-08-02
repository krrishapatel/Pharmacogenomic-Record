from pathlib import Path

from pgxrecord.ingest.raw import RawCall
from pgxrecord.ingest.vcf import build_vcf
from pgxrecord.positions import ReferencePosition

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
    from pgxrecord.ingest.vcf import translate_genotype

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


def test_gene_with_no_rsid_positions_is_always_fully_uncovered(tmp_path):
    """CYP2D6 relies on positions with no rsID, so an array can never cover it.

    This is the single most important coverage case: CYP2D6 is among the most
    clinically significant PGx genes and consumer arrays cannot resolve it.
    """
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="GG")]

    report = build_vcf(calls, REF, out)

    assert "CYP2D6" in report.genes_fully_uncovered


def test_untranslatable_genotype_counts_as_uncovered(tmp_path):
    """A call whose alleles don't match ref/alt yields no data, not a ref call."""
    out = tmp_path / "out.vcf"
    calls = [RawCall(rsid="rs1", chrom="1", pos=100, genotype="AA")]

    report = build_vcf(calls, REF, out)

    assert "rs1" not in report.covered_rsids
    assert "rs1" in report.uncovered_rsids
    assert "rs1" not in out.read_text()
