from pathlib import Path

import pytest

from pharmacogenomic_record import POSITIONS_FILENAME
from dataclasses import FrozenInstanceError

from pharmacogenomic_record.positions import (
    ReferencePosition,
    genes_covered,
    index_by_rsid,
    load_positions,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
POSITIONS = REPO_ROOT / "src" / "pharmacogenomic_record" / "data" / POSITIONS_FILENAME


@pytest.fixture(scope="module")
def positions():
    return load_positions(POSITIONS)


def test_load_positions_parses_every_data_row(positions):
    assert len(positions) == 1226


def test_parsed_fields_match_the_file(positions):
    first = positions[0]
    assert first.chrom == "chr1"
    assert first.pos == 97078987
    assert first.rsid == "rs114096998"
    assert first.ref == "G"
    assert first.alt == ("T",)
    assert first.gene == "DPYD"


def test_multi_allelic_alt_is_split(positions):
    """57 of 1226 positions are multi-allelic; ALT must be split on commas.

    Without this, a single test on a single-allelic row lets a broken
    implementation (alt=(raw,)) pass, and downstream genotype matching in
    the ingest step would silently fail on every multi-allelic position.
    """
    by_rsid = {p.rsid: p for p in positions if p.rsid}
    assert by_rsid["rs3064744"].ref == "CAT"
    assert by_rsid["rs3064744"].alt == ("C", "CATAT", "CATATAT")
    assert len([p for p in positions if len(p.alt) > 1]) == 57


def test_reference_position_is_hashable(positions):
    """A tuple alt keeps positions usable in sets and as dict keys."""
    assert len({p for p in positions}) == 1226


def test_positions_without_rsid_get_none(positions):
    """208 positions in 3.4.0 have '.' as ID. These are never joinable."""
    without = [p for p in positions if p.rsid is None]
    assert len(without) == 208


def test_index_by_rsid_skips_unjoinable_positions(positions):
    index = index_by_rsid(positions)
    assert len(index) == 1018
    assert "rs114096998" in index
    assert index["rs114096998"].gene == "DPYD"
    assert None not in index
    assert "." not in index


def test_position_without_gene_tag_is_kept_with_gene_none(positions):
    """rs12777823 has INFO 'POI' and no PX= tag -- 1225 of 1226 have a gene.

    It must parse, not raise, and must not be counted as a gene.
    """
    by_rsid = {p.rsid: p for p in positions if p.rsid}
    assert by_rsid["rs12777823"].gene is None
    assert len([p for p in positions if p.gene is not None]) == 1225


def test_genes_covered(positions):
    genes = genes_covered(positions)
    assert len(genes) == 22
    assert {"CYP2C19", "CYP2D6", "DPYD", "TPMT", "SLCO1B1"} <= genes


def test_reference_position_is_immutable():
    p = ReferencePosition(
        chrom="chr1", pos=1, rsid="rs1", ref="A", alt=("G",), gene="DPYD"
    )
    with pytest.raises(FrozenInstanceError):
        p.pos = 2


def test_malformed_line_error_names_the_file_and_line_number(tmp_path):
    """A row that cannot be split into 8 columns must locate itself.

    A bare "not enough values to unpack" names neither file nor line, so a user
    with a 1200-line reference table has nothing to go on. The message must
    carry the path and the 1-based line number, matching ingest/raw.py's
    `{path}:{number}:` style.
    """
    bad = tmp_path / "positions.vcf"
    # A valid header line, one good data row, then a truncated data row on line 3.
    bad.write_text(
        "##fileformat=VCFv4.2\n"
        "chr1\t1\trs1\tA\tG\t.\tPASS\tPX=DPYD\n"
        "chr1\t2\trs2\n"
    )

    with pytest.raises(ValueError) as excinfo:
        load_positions(bad)

    message = str(excinfo.value)
    assert str(bad) in message
    # The truncated row is line 3 (1-based, counting the header).
    assert ":3:" in message
