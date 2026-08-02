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
