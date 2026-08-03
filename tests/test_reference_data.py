from pathlib import Path

from pharmacogenomic_record import POSITIONS_FILENAME

REPO_ROOT = Path(__file__).resolve().parents[1]
POSITIONS = REPO_ROOT / "src" / "pharmacogenomic_record" / "data" / POSITIONS_FILENAME


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
    """Confirms the build, which is why we join on rsID rather than position.

    Asserts on the contig header so the exact patch level is pinned, not just
    the substring "GRCh38" appearing somewhere in the file.
    """
    assert '##contig=<ID=chr1,assembly=GRCh38.p14' in POSITIONS.read_text()
