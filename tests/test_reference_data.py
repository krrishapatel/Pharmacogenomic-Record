from pathlib import Path

from pgxrecord import POSITIONS_FILENAME

POSITIONS = Path("data") / POSITIONS_FILENAME


def test_positions_file_exists():
    assert POSITIONS.is_file()


def test_positions_file_has_expected_shape():
    """The reference file must match the pinned PharmCAT version exactly.

    These counts are from pharmcat_positions_3.4.0.vcf. If they change, the
    pinned version changed and every stored record's guideline_version stamp
    is now suspect.
    """
    data_lines = [
        line
        for line in POSITIONS.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(data_lines) == 1226

    with_rsid = [line for line in data_lines if line.split("\t")[2].startswith("rs")]
    assert len(with_rsid) == 1018

    genes = {
        field[3:]
        for line in data_lines
        for field in line.split("\t")[7].split(";")
        if field.startswith("PX=")
    }
    assert "CYP2C19" in genes
    assert "DPYD" in genes
    assert len(genes) == 22


def test_positions_file_is_grch38():
    """Confirms the build, which is why we join on rsID rather than position."""
    assert "GRCh38" in POSITIONS.read_text(errors="ignore")[:20000]
