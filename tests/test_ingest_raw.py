from pathlib import Path

import pytest

from pgxrecord.ingest.raw import RawCall, UnsupportedRawFile, parse_23andme

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_parses_genotype_rows():
    calls = parse_23andme(FIXTURES / "23andme_valid.txt")
    by_rsid = {c.rsid: c for c in calls}
    assert by_rsid["rs4244285"] == RawCall(
        rsid="rs4244285", chrom="10", pos=96541616, genotype="AG"
    )


def test_skips_internal_ids_and_nocalls():
    """Internal 'i' IDs are unjoinable; '--' means the array failed to call."""
    calls = parse_23andme(FIXTURES / "23andme_valid.txt")
    rsids = {c.rsid for c in calls}
    assert "i5000123" not in rsids
    assert "rs9999999999" not in rsids
    assert len(calls) == 4


def test_rejects_non_build37_file():
    """Build 38 raw files would silently break the rsID join assumption."""
    with pytest.raises(UnsupportedRawFile, match="build 37"):
        parse_23andme(FIXTURES / "23andme_build38.txt")


def test_rejects_file_without_recognizable_header():
    """Never guess the vendor or build. Reject instead."""
    with pytest.raises(UnsupportedRawFile, match="header"):
        parse_23andme(FIXTURES / "23andme_no_header.txt")


def test_rejects_build38_file_that_mentions_build_37_in_prose():
    """The build must be read positionally, not by substring search.

    A substring check for "build 37" passes this file, and every coordinate
    downstream would then be GRCh38 data treated as GRCh37 -- silent, wrong
    star-allele calls. This is the failure the rsID join exists to avoid.
    """
    with pytest.raises(UnsupportedRawFile, match="more than one reference"):
        parse_23andme(FIXTURES / "23andme_build38_mentions_37.txt")


def test_rejects_unambiguous_build38_file():
    """A header declaring only build 38 is rejected on the build number."""
    with pytest.raises(UnsupportedRawFile, match="declares build 38"):
        parse_23andme(FIXTURES / "23andme_build38.txt")


def test_rejects_wrong_column_count_instead_of_returning_nothing():
    """A wrong-format file must fail loudly, not parse to zero calls.

    Silently skipping rows of the wrong arity means an AncestryDNA or
    MyHeritage export with a compatible-looking header returns [] and is
    indistinguishable from a clean file containing nothing relevant.
    """
    with pytest.raises(UnsupportedRawFile, match="expected 4 tab-separated"):
        parse_23andme(FIXTURES / "23andme_five_columns.txt")


def test_rejects_non_numeric_position():
    """int(pos) must not escape the module's UnsupportedRawFile contract."""
    bad = FIXTURES / "23andme_bad_position.txt"
    bad.write_text(
        "# 23andMe\n"
        "# We are using reference human assembly build 37.\n"
        "rs114096998\t1\tNOTANUMBER\tGG\n"
    )
    try:
        with pytest.raises(UnsupportedRawFile, match="is not a number"):
            parse_23andme(bad)
    finally:
        bad.unlink()


def test_drops_hemizygous_indel_and_single_dash_genotypes():
    """Only diploid nucleotide genotypes survive.

    A 1-character genotype has no valid diploid VCF GT, a single dash is a
    hemizygous no-call, and D/I indel codes are not nucleotide alleles that
    can be matched against reference ref/alt bases.
    """
    calls = parse_23andme(FIXTURES / "23andme_edge_genotypes.txt")

    assert {c.rsid for c in calls} == {"rs114096998", "rs1801268"}
    assert all(len(c.genotype) == 2 for c in calls)
    assert all(not set(c.genotype) & {"D", "I"} for c in calls)


def test_rejects_file_that_yields_no_usable_calls():
    """Parsing to an empty list is a format failure, not a clean result."""
    empty = FIXTURES / "23andme_all_nocalls.txt"
    empty.write_text(
        "# 23andMe\n"
        "# We are using reference human assembly build 37.\n"
        "rs114096998\t1\t97544543\t--\n"
    )
    try:
        with pytest.raises(UnsupportedRawFile, match="no usable genotype"):
            parse_23andme(empty)
    finally:
        empty.unlink()
