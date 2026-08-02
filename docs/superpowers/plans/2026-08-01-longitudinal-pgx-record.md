# Longitudinal PGx Record Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ingest a 23andMe raw genotype file, call pharmacogenomic star alleles via PharmCAT, and store the result as an immutable version-stamped record that can be re-queried per drug and re-evaluated when guidelines change.

**Architecture:** Four isolated Python packages — `ingest` (raw text → VCF, joined on rsID), `call` (thin subprocess wrapper around the pinned PharmCAT Docker image), `record` (append-only SQLite store), `evaluate` (drug query + guideline-bump diffing). Data flows one direction; genotype is called once and queried many times.

**Tech Stack:** Python 3.11+, pytest, SQLite (stdlib `sqlite3`), Docker (`pgkb/pharmcat:3.4.0`). No web framework in v1 — library plus CLI only.

## Global Constraints

- **PharmCAT is pinned to `pgkb/pharmcat:3.4.0`.** Never `latest`. Never reimplement allele calling.
- **Reference positions come from `pharmcat_positions_3.4.0.vcf`** (1,226 positions, 22 genes, GRCh38.p14) — downloaded once into `data/`, committed to the repo.
- **Join genotype to reference positions on rsID, never on coordinates.** 23andMe raw data is GRCh37; PharmCAT positions are GRCh38. rsID join avoids liftover entirely.
- **Guideline content is referenced by identifier and URL, never stored as prose.**
- **Three coverage states, never collapsed:** `called`, `not_covered`, `indeterminate`.
- **Three query outcomes, never collapsed:** `guidance_found`, `no_guidance_for_pair`, `cannot_assess`.
- **`cannot_assess` must never render as reassurance.** No code path may return an empty/falsy "nothing found" for an unassessable gene.
- **Records are append-only.** No `UPDATE`, no `DELETE` on the records table, ever.
- **No clinical claims in any output string.** No dose figures, no "you should" phrasing.
- **Docker is NOT installed on the development machine.** Tasks 1-8 are entirely
  Docker-free and must stay that way. Only Task 9's manual end-to-end ingest
  needs it, and that step is documented as unverified rather than skipped.
- **Genotype data is never committed.** `data/.gitignore` denies everything
  except the two committed reference tables; PharmCAT scratch output goes to
  `work/`, which is git-ignored.

---

### Task 1: Project scaffold and reference position data

**Files:**
- Create: `pyproject.toml`
- Create: `src/pgxrecord/__init__.py`
- Create: `data/.gitignore`
- Create: `scripts/fetch_positions.sh`
- Test: `tests/test_reference_data.py`

**Interfaces:**
- Consumes: nothing (first task)
- Produces: package `pgxrecord` importable; `data/pharmcat_positions_3.4.0.vcf` present on disk; constant `PHARMCAT_VERSION = "3.4.0"` in `src/pgxrecord/__init__.py`

- [ ] **Step 1: Create the package scaffold**

`pyproject.toml`:

```toml
[project]
name = "pgxrecord"
version = "0.1.0"
description = "Longitudinal pharmacogenomic record over PharmCAT"
requires-python = ">=3.11"
dependencies = []

[project.optional-dependencies]
dev = ["pytest>=8.0"]

[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

`src/pgxrecord/__init__.py`:

```python
"""Longitudinal pharmacogenomic record built on PharmCAT."""

PHARMCAT_VERSION = "3.4.0"
PHARMCAT_IMAGE = f"pgkb/pharmcat:{PHARMCAT_VERSION}"
POSITIONS_FILENAME = f"pharmcat_positions_{PHARMCAT_VERSION}.vcf"
```

`data/.gitignore`. The reference position table and the gene-drug pair table
are committed; everything else that lands here is derived or personal
genotype data and must never be committed:

```
# Ignore everything by default -- genotype data must never be committed.
*
# ...except this file and the two committed reference tables.
!.gitignore
!pharmcat_positions_*.vcf
!gene_drug_pairs.json
```

- [ ] **Step 2: Write the fetch script**

`scripts/fetch_positions.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail
VERSION="3.4.0"
EXPECTED_BYTES=64934
DEST="data/pharmcat_positions_${VERSION}.vcf"
URL="https://github.com/PharmGKB/PharmCAT/releases/download/v${VERSION}/pharmcat_positions_${VERSION}.vcf"
mkdir -p data

# Download to a temp file and validate before replacing the pinned reference.
# Without -f, curl writes the HTTP error body to the output file and exits 0,
# which would silently overwrite the good file with "Not Found".
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
curl -fsSL -o "$TMP" "$URL"

ACTUAL_BYTES="$(wc -c < "$TMP" | tr -d " ")"
if [ "$ACTUAL_BYTES" -ne "$EXPECTED_BYTES" ]; then
    echo "refusing to install: expected ${EXPECTED_BYTES} bytes, got ${ACTUAL_BYTES}" >&2
    exit 1
fi

mv "$TMP" "$DEST"
trap - EXIT
echo "wrote $DEST (${ACTUAL_BYTES} bytes)"
```

Run:

```bash
chmod +x scripts/fetch_positions.sh && ./scripts/fetch_positions.sh
```

Expected: `wrote data/pharmcat_positions_3.4.0.vcf (64934 bytes)`

- [ ] **Step 3: Write the failing test**

`tests/test_reference_data.py`:

```python
from pathlib import Path

from pgxrecord import POSITIONS_FILENAME

REPO_ROOT = Path(__file__).resolve().parents[1]
POSITIONS = REPO_ROOT / "data" / POSITIONS_FILENAME


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
```

- [ ] **Step 4: Run tests to verify they fail, then pass**

Run: `pip install -e ".[dev]" && pytest tests/test_reference_data.py -v`

Expected: PASS (3 tests). If `test_positions_file_exists` fails, Step 2 was not run.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/pgxrecord/__init__.py data/.gitignore scripts/fetch_positions.sh tests/test_reference_data.py data/pharmcat_positions_3.4.0.vcf
git commit -m "feat: scaffold package and pin PharmCAT 3.4.0 reference positions"
```

---

### Task 2: Parse the reference position table

**Files:**
- Create: `src/pgxrecord/positions.py`
- Test: `tests/test_positions.py`

**Interfaces:**
- Consumes: `data/pharmcat_positions_3.4.0.vcf` from Task 1
- Produces:
  - `class ReferencePosition` — frozen dataclass with fields `chrom: str`, `pos: int`, `rsid: str | None`, `ref: str`, `alt: tuple[str, ...]`, `gene: str | None`
  - `def load_positions(path: Path) -> list[ReferencePosition]`
  - `def index_by_rsid(positions: list[ReferencePosition]) -> dict[str, ReferencePosition]`
  - `def genes_covered(positions: list[ReferencePosition]) -> set[str]`

- [ ] **Step 1: Write the failing test**

`tests/test_positions.py`:

```python
from pathlib import Path

import pytest

from pgxrecord import POSITIONS_FILENAME
from dataclasses import FrozenInstanceError

from pgxrecord.positions import (
    ReferencePosition,
    genes_covered,
    index_by_rsid,
    load_positions,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
POSITIONS = REPO_ROOT / "data" / POSITIONS_FILENAME


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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_positions.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.positions'`

- [ ] **Step 3: Write the implementation**

`src/pgxrecord/positions.py`:

```python
"""Parse the PharmCAT reference position table.

The reference file is the authoritative list of positions PharmCAT can call,
shipped with the pinned PharmCAT release. Positions whose ID is '.' have no
rsID and therefore cannot be joined against consumer array data at all --
they always become not_covered.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReferencePosition:
    """A single position PharmCAT knows how to interpret."""

    chrom: str
    pos: int
    rsid: str | None
    ref: str
    alt: tuple[str, ...]
    gene: str | None


def _parse_gene(info: str) -> str | None:
    """Extract the PX= gene tag, or None when the position has no gene.

    Exactly one position in 3.4.0 (rs12777823, chr10:94645745) carries INFO
    'POI' -- a position of interest with no gene assignment. It is a real,
    joinable position, so we keep it and leave gene as None rather than
    rejecting the file.
    """
    for field in info.split(";"):
        if field.startswith("PX="):
            return field[3:]
    return None


def load_positions(path: Path) -> list[ReferencePosition]:
    """Parse every data row of a PharmCAT positions VCF."""
    positions: list[ReferencePosition] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        chrom, pos, rsid, ref, alt, _qual, _filter, info = line.split("\t")[:8]
        positions.append(
            ReferencePosition(
                chrom=chrom,
                pos=int(pos),
                rsid=rsid if rsid.startswith("rs") else None,
                ref=ref,
                alt=tuple(alt.split(",")),
                gene=_parse_gene(info),
            )
        )
    return positions


def index_by_rsid(
    positions: list[ReferencePosition],
) -> dict[str, ReferencePosition]:
    """Index joinable positions by rsID.

    Positions without an rsID are omitted. We join on rsID rather than
    coordinates because consumer arrays report GRCh37 while this file is
    GRCh38; rsID avoids a liftover step and its attendant errors.
    """
    return {p.rsid: p for p in positions if p.rsid is not None}


def genes_covered(positions: list[ReferencePosition]) -> set[str]:
    """Return every gene appearing in the reference table.

    Positions with no gene assignment (INFO 'POI') are excluded.
    """
    return {p.gene for p in positions if p.gene is not None}
```

Note: `alt` is a `tuple`, not a `list`. `frozen=True` blocks attribute
rebinding but does NOT stop `position.alt.append(...)` on a contained list,
and a list field also makes `ReferencePosition` unhashable — so it could
never be put in a `set` or used as a dict key by a later task. A tuple fixes
both.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_positions.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pgxrecord/positions.py tests/test_positions.py
git commit -m "feat: parse PharmCAT reference position table, indexed by rsID"
```

---

### Task 3: Parse 23andMe raw genotype files

**Files:**
- Create: `src/pgxrecord/ingest/__init__.py`
- Create: `src/pgxrecord/ingest/raw.py`
- Create: `tests/fixtures/23andme_valid.txt`
- Create: `tests/fixtures/23andme_build38.txt`
- Create: `tests/fixtures/23andme_no_header.txt`
- Test: `tests/test_ingest_raw.py`

**Interfaces:**
- Consumes: nothing from prior tasks
- Produces:
  - `class RawCall` — frozen dataclass, fields `rsid: str`, `chrom: str`, `pos: int`, `genotype: str`
  - `class UnsupportedRawFile(Exception)`
  - `def parse_23andme(path: Path) -> list[RawCall]` — raises `UnsupportedRawFile` when the build is not 37 or the header is absent

- [ ] **Step 1: Create the fixtures**

`tests/fixtures/23andme_valid.txt` (tab-separated; `rs1801268` and `rs114096998` are real DPYD positions from the reference table, `rs4244285` is CYP2C19*2):

```
# This data file generated by 23andMe at: Sat Jan  1 00:00:00 2022
# This file contains raw genotype data, including data that is not used in 23andMe reports.
# We are using reference human assembly build 37 (also known as Annotation Release 104).
# rsid	chromosome	position	genotype
rs114096998	1	97544543	GG
rs1801268	1	97544627	CC
rs4244285	10	96541616	AG
rs12248560	10	96521657	CC
i5000123	1	12345	AA
rs9999999999	1	99999	--
```

`tests/fixtures/23andme_build38.txt`:

```
# This data file generated by 23andMe at: Sat Jan  1 00:00:00 2022
# We are using reference human assembly build 38.
# rsid	chromosome	position	genotype
rs114096998	1	97078987	GG
```

`tests/fixtures/23andme_no_header.txt`:

```
rs114096998	1	97544543	GG
rs1801268	1	97544627	CC
```

- [ ] **Step 2: Write the failing test**

`tests/test_ingest_raw.py`:

```python
from pathlib import Path

import pytest

from pgxrecord.ingest.raw import RawCall, UnsupportedRawFile, parse_23andme

FIXTURES = Path("tests/fixtures")


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
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_ingest_raw.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.ingest'`

- [ ] **Step 4: Write the implementation**

`src/pgxrecord/ingest/__init__.py`:

```python
"""Convert consumer raw genotype files into PharmCAT-ready VCF."""
```

`src/pgxrecord/ingest/raw.py`:

```python
"""Parse 23andMe raw genotype exports.

Consumer array files are messy in specific known ways, and this module is
where that mess is contained. Two rules:

1. Never guess the vendor or genome build. An unrecognized header is a
   rejection, not a default.
2. Drop rows that cannot be joined -- internal 'i' identifiers and no-calls
   ('--'). Dropping them here means downstream code sees only usable calls,
   and the positions they would have covered surface as not_covered.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

NO_CALL = "--"


class UnsupportedRawFile(Exception):
    """The raw file is not a format we can convert safely."""


@dataclass(frozen=True)
class RawCall:
    """One genotype call from a consumer array."""

    rsid: str
    chrom: str
    pos: int
    genotype: str


def _validate_header(header_text: str) -> None:
    lowered = header_text.lower()
    if "23andme" not in lowered:
        raise UnsupportedRawFile(
            "no recognizable 23andMe header found; refusing to guess the "
            "vendor or genome build"
        )
    if "build 37" not in lowered:
        raise UnsupportedRawFile(
            "expected reference assembly build 37; this file declares a "
            "different build, which would invalidate the rsID join"
        )


def parse_23andme(path: Path) -> list[RawCall]:
    """Parse a 23andMe raw export into joinable genotype calls."""
    lines = path.read_text().splitlines()
    header_text = "\n".join(line for line in lines if line.startswith("#"))
    _validate_header(header_text)

    calls: list[RawCall] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 4:
            continue
        rsid, chrom, pos, genotype = fields
        if not rsid.startswith("rs"):
            continue
        if genotype == NO_CALL:
            continue
        calls.append(
            RawCall(rsid=rsid, chrom=chrom, pos=int(pos), genotype=genotype)
        )
    return calls
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_ingest_raw.py -v`
Expected: PASS (4 tests)

- [ ] **Step 6: Commit**

```bash
git add src/pgxrecord/ingest tests/fixtures tests/test_ingest_raw.py
git commit -m "feat: parse 23andMe raw exports, rejecting unknown vendor and build"
```

---

### Task 4: Build VCF from raw calls, with an explicit coverage report

**Files:**
- Create: `src/pgxrecord/ingest/vcf.py`
- Test: `tests/test_ingest_vcf.py`

**Interfaces:**
- Consumes: `RawCall` (Task 3), `ReferencePosition` / `load_positions` / `index_by_rsid` (Task 2)
- Produces:
  - `class CoverageReport` — frozen dataclass, fields `covered_rsids: set[str]`, `uncovered_rsids: set[str]`, `genes_fully_uncovered: set[str]`, `genes_partially_covered: set[str]`
  - `def build_vcf(calls: list[RawCall], positions: list[ReferencePosition], out_path: Path) -> CoverageReport`

This is the task that makes the `cannot_assess` invariant possible. If coverage is not computed here, nothing downstream can distinguish "no interaction" from "no data."

- [ ] **Step 1: Write the failing test**

`tests/test_ingest_vcf.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ingest_vcf.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.ingest.vcf'`

- [ ] **Step 3: Write the implementation**

`src/pgxrecord/ingest/vcf.py`:

```python
"""Emit a PharmCAT-ready VCF from consumer array calls.

Two things make this module load-bearing for correctness:

1. Coordinates always come from the reference table (GRCh38), never from the
   raw file (GRCh37). The rsID is the join key; the raw position is discarded.
2. The CoverageReport is what lets downstream code answer "we do not know"
   instead of "no interaction found". A gene whose positions are absent from
   the array carries no information, and that must stay visible.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pgxrecord.ingest.raw import RawCall
from pgxrecord.positions import ReferencePosition, index_by_rsid

_VCF_META = """##fileformat=VCFv4.2
##source=pgxrecord
##reference=GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
_VCF_COLUMNS = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Order chromosomes naturally: chr1, chr2, ... chr10, ... chrX, chrY.

    Plain string sort puts chr10 before chr2, which produces an out-of-order
    VCF. Numeric contigs sort by value; X/Y/M sort after them.
    """
    name = chrom.removeprefix("chr")
    return (int(name), "") if name.isdigit() else (10**6, name)


def _contig_lines(positions: list[ReferencePosition]) -> str:
    """Declare every contig we emit, in sorted order.

    VCF consumers may reject or misparse records whose contig was never
    declared in the header.
    """
    chroms = sorted({p.chrom for p in positions}, key=_chrom_sort_key)
    return "".join(
        f'##contig=<ID={chrom},assembly=GRCh38.p14,species="Homo sapiens">\n'
        for chrom in chroms
    )


@dataclass(frozen=True)
class CoverageReport:
    """Which reference positions the array actually informed."""

    covered_rsids: set[str]
    uncovered_rsids: set[str]
    genes_fully_uncovered: set[str]
    genes_partially_covered: set[str]


def translate_genotype(genotype: str, ref: ReferencePosition) -> str | None:
    """Convert raw allele letters to a VCF numeric genotype.

    Returns None when any allele is neither the reference nor a known
    alternate, which means the call tells us nothing about this position.
    """
    if len(genotype) != 2:
        return None
    alleles = [ref.ref, *ref.alt]
    try:
        indices = sorted(alleles.index(base) for base in genotype)
    except ValueError:
        return None
    return f"{indices[0]}/{indices[1]}"


def build_vcf(
    calls: list[RawCall],
    positions: list[ReferencePosition],
    out_path: Path,
) -> CoverageReport:
    """Write a VCF for every reference position the array covers."""
    by_rsid = index_by_rsid(positions)
    calls_by_rsid = {c.rsid: c for c in calls}

    rows: list[str] = []
    covered: set[str] = set()

    for rsid, ref in by_rsid.items():
        call = calls_by_rsid.get(rsid)
        if call is None:
            continue
        gt = translate_genotype(call.genotype, ref)
        if gt is None:
            continue
        covered.add(rsid)
        rows.append(
            f"{ref.chrom}\t{ref.pos}\t{ref.rsid}\t{ref.ref}\t"
            f"{','.join(ref.alt)}\t.\tPASS\t.\tGT\t{gt}"
        )

    rows.sort(
        key=lambda row: (
            _chrom_sort_key(row.split("\t")[0]),
            int(row.split("\t")[1]),
        )
    )
    out_path.write_text(
        _VCF_META
        + _contig_lines(positions)
        + _VCF_COLUMNS
        + "".join(f"{row}\n" for row in rows)
    )

    all_joinable = set(by_rsid)
    # Positions with no PX= gene tag (INFO 'POI') contribute no gene.
    genes_covered_partly = {
        by_rsid[r].gene for r in covered if by_rsid[r].gene is not None
    }
    all_genes = {p.gene for p in positions if p.gene is not None}

    return CoverageReport(
        covered_rsids=covered,
        uncovered_rsids=all_joinable - covered,
        genes_fully_uncovered=all_genes - genes_covered_partly,
        genes_partially_covered=genes_covered_partly,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ingest_vcf.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pgxrecord/ingest/vcf.py tests/test_ingest_vcf.py
git commit -m "feat: build VCF from raw calls with explicit coverage reporting"
```

---

### Task 5: PharmCAT subprocess wrapper

**Files:**
- Create: `src/pgxrecord/caller.py`
- Create: `tests/fixtures/pharmcat_phenotype_sample.json`
- Test: `tests/test_caller.py`

**Interfaces:**
- Consumes: `PHARMCAT_IMAGE` (Task 1)
- Produces:
  - `class GeneCall` — frozen dataclass, fields `gene: str`, `diplotype: str | None`, `phenotype: str | None`, `coverage: str` (one of `"called"`, `"not_covered"`, `"indeterminate"`)
  - `class PharmcatError(Exception)`
  - `def run_pharmcat(vcf_path: Path, workdir: Path) -> Path` — returns path to `*.phenotype.json`; raises `PharmcatError`
  - `def parse_phenotype_json(path: Path, uncovered_genes: set[str]) -> list[GeneCall]`

We never implement calling logic. This wrapper only invokes the pinned image and parses its output.

- [ ] **Step 1: Create the fixture**

`tests/fixtures/pharmcat_phenotype_sample.json` — a trimmed but structurally faithful PharmCAT phenotype output:

```json
{
  "phenotypes": {
    "CYP2C19": {
      "gene": "CYP2C19",
      "diplotypes": [
        {
          "allele1": {"name": "*1"},
          "allele2": {"name": "*2"},
          "label": "*1/*2",
          "phenotypes": ["Intermediate Metabolizer"]
        }
      ]
    },
    "DPYD": {
      "gene": "DPYD",
      "diplotypes": [
        {
          "allele1": {"name": "Reference"},
          "allele2": {"name": "Reference"},
          "label": "Reference/Reference",
          "phenotypes": ["Normal Metabolizer"]
        }
      ]
    },
    "TPMT": {
      "gene": "TPMT",
      "diplotypes": [
        {
          "allele1": null,
          "allele2": null,
          "label": "Unknown/Unknown",
          "phenotypes": ["Indeterminate"]
        }
      ]
    }
  }
}
```

- [ ] **Step 2: Write the failing test**

`tests/test_caller.py`:

```python
from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall, PharmcatError, parse_phenotype_json, run_pharmcat

FIXTURES = Path("tests/fixtures")
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


def test_malformed_json_raises():
    bad = FIXTURES / "pharmcat_malformed.json"
    bad.write_text("{not json")
    try:
        with pytest.raises(PharmcatError, match="parse"):
            parse_phenotype_json(bad, uncovered_genes=set())
    finally:
        bad.unlink()


def test_run_pharmcat_raises_when_docker_missing(tmp_path, monkeypatch):
    """No record may be written when the caller cannot run."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    monkeypatch.setenv("PATH", str(tmp_path))  # hide docker

    with pytest.raises(PharmcatError):
        run_pharmcat(vcf, tmp_path)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_caller.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.caller'`

- [ ] **Step 4: Write the implementation**

`src/pgxrecord/caller.py`:

```python
"""Thin wrapper around the pinned PharmCAT Docker image.

We do not implement star-allele calling and never will -- PharmCAT is
maintained by the CPIC/PharmGKB group and is the authority. This module only
invokes it and translates its output into our types.

Two rules:

1. The image tag is pinned. A guideline or tool update must be a deliberate,
   reviewed act, because every stored record is stamped with the version that
   produced it.
2. If PharmCAT fails or emits anything unparseable, raise. A partial record
   is worse than no record.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from pgxrecord import PHARMCAT_IMAGE

CALLED = "called"
NOT_COVERED = "not_covered"
INDETERMINATE = "indeterminate"

_INDETERMINATE_LABELS = {"indeterminate", "unknown", "n/a", "no result"}
_TIMEOUT_SECONDS = 900


class PharmcatError(Exception):
    """PharmCAT could not be run, or its output could not be trusted."""


@dataclass(frozen=True)
class GeneCall:
    """One gene's result, with explicit coverage state."""

    gene: str
    diplotype: str | None
    phenotype: str | None
    coverage: str


def run_pharmcat(vcf_path: Path, workdir: Path) -> Path:
    """Run the pinned PharmCAT image over a VCF and return the phenotype JSON."""
    if shutil.which("docker") is None:
        raise PharmcatError("docker not found on PATH; cannot run PharmCAT")

    workdir = workdir.resolve()
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{workdir}:/pharmcat/data",
        PHARMCAT_IMAGE,
        "pharmcat_pipeline", f"/pharmcat/data/{vcf_path.name}",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=_TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired as err:
        raise PharmcatError(f"PharmCAT timed out after {_TIMEOUT_SECONDS}s") from err
    except OSError as err:
        raise PharmcatError(f"could not execute docker: {err}") from err

    if result.returncode != 0:
        raise PharmcatError(
            f"PharmCAT exited {result.returncode}\nstderr:\n{result.stderr}"
        )

    phenotype = workdir / f"{vcf_path.stem}.phenotype.json"
    if not phenotype.is_file():
        raise PharmcatError(
            f"PharmCAT produced no phenotype JSON at {phenotype}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return phenotype


def _classify(label: str | None, phenotype: str | None) -> str:
    candidates = [value.lower() for value in (label, phenotype) if value]
    if not candidates:
        return INDETERMINATE
    if any(
        marker in value
        for value in candidates
        for marker in _INDETERMINATE_LABELS
    ):
        return INDETERMINATE
    return CALLED


def parse_phenotype_json(path: Path, uncovered_genes: set[str]) -> list[GeneCall]:
    """Translate PharmCAT phenotype output into GeneCall records.

    uncovered_genes always wins. If the array never informed a gene, no
    PharmCAT output about it can be treated as a call.
    """
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as err:
        raise PharmcatError(f"could not parse PharmCAT output {path}: {err}") from err

    calls: list[GeneCall] = []
    seen: set[str] = set()

    for gene, entry in (payload.get("phenotypes") or {}).items():
        seen.add(gene)
        if gene in uncovered_genes:
            calls.append(GeneCall(gene, None, None, NOT_COVERED))
            continue

        diplotypes = entry.get("diplotypes") or [{}]
        first = diplotypes[0]
        label = first.get("label")
        phenotypes = first.get("phenotypes") or []
        phenotype = phenotypes[0] if phenotypes else None
        coverage = _classify(label, phenotype)

        calls.append(
            GeneCall(
                gene=gene,
                diplotype=label if coverage == CALLED else None,
                phenotype=phenotype if coverage == CALLED else None,
                coverage=coverage,
            )
        )

    for gene in sorted(uncovered_genes - seen):
        calls.append(GeneCall(gene, None, None, NOT_COVERED))

    return calls
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_caller.py -v`
Expected: PASS (6 tests)

- [ ] **Step 6: Commit**

```bash
git add src/pgxrecord/caller.py tests/fixtures/pharmcat_phenotype_sample.json tests/test_caller.py
git commit -m "feat: wrap pinned PharmCAT image and parse phenotype output"
```

---

### Task 6: Append-only record store

**Files:**
- Create: `src/pgxrecord/store.py`
- Test: `tests/test_store.py`

**Interfaces:**
- Consumes: `GeneCall` (Task 5), `PHARMCAT_VERSION` (Task 1)
- Produces:
  - `class RecordStore` with methods:
    - `__init__(self, db_path: Path)`
    - `append(self, subject_id: str, calls: list[GeneCall], guideline_version: str) -> int` — returns `record_id`
    - `latest(self, subject_id: str) -> list[GeneCall]`
    - `history(self, subject_id: str) -> list[int]` — record_ids oldest to newest
    - `record_versions(self, record_id: int) -> tuple[str, str]` — `(pharmcat_version, guideline_version)`
    - `subjects_with_gene(self, gene: str) -> list[str]`

- [ ] **Step 1: Write the failing test**

`tests/test_store.py`:

```python
import sqlite3

import pytest

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import GeneCall
from pgxrecord.store import RecordStore

CALLS = [
    GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
    GeneCall("CYP2D6", None, None, "not_covered"),
]


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


def test_append_and_read_back(store):
    store.append("subject-1", CALLS, guideline_version="cpic-2026-07")
    calls = store.latest("subject-1")

    by_gene = {c.gene: c for c in calls}
    assert by_gene["CYP2C19"].phenotype == "Intermediate Metabolizer"
    assert by_gene["CYP2D6"].coverage == "not_covered"


def test_versions_are_stamped(store):
    record_id = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    pharmcat_version, guideline_version = store.record_versions(record_id)

    assert pharmcat_version == PHARMCAT_VERSION
    assert guideline_version == "cpic-2026-07"


def test_reingest_appends_and_preserves_history(store):
    first = store.append("s1", CALLS, guideline_version="cpic-2026-07")
    updated = [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called")]
    second = store.append("s1", updated, guideline_version="cpic-2026-08")

    assert store.history("s1") == [first, second]
    assert store.latest("s1")[0].phenotype == "Normal Metabolizer"
    # The original call is still retrievable -- history is not destroyed.
    assert store.record_versions(first)[1] == "cpic-2026-07"


def test_updates_are_rejected_at_the_database_level(store):
    """Immutability is enforced by the schema, not by convention."""
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("UPDATE gene_calls SET phenotype = 'tampered'")
        conn.commit()
    conn.close()


def test_deletes_are_rejected_at_the_database_level(store):
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    conn = sqlite3.connect(store.db_path)
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("DELETE FROM gene_calls")
        conn.commit()
    conn.close()


def test_subjects_with_gene_finds_affected_records(store):
    store.append("s1", CALLS, guideline_version="cpic-2026-07")
    store.append("s2", [GeneCall("DPYD", "Ref/Ref", "Normal", "called")],
                 guideline_version="cpic-2026-07")

    assert store.subjects_with_gene("CYP2C19") == ["s1"]
    assert store.subjects_with_gene("DPYD") == ["s2"]


def test_unknown_subject_returns_empty_history(store):
    assert store.history("nobody") == []
    assert store.latest("nobody") == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_store.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.store'`

- [ ] **Step 3: Write the implementation**

`src/pgxrecord/store.py`:

```python
"""Append-only pharmacogenomic record store.

Immutability is a correctness requirement, not a preference. A phenotype call
only means something relative to the tool and guideline versions that produced
it. Overwriting a call destroys the ability to explain why the system once
said something different -- which is exactly the question that matters when
guidance is revised.

SQLite triggers enforce this so that a future careless UPDATE fails loudly
rather than silently rewriting clinical history.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from pgxrecord import PHARMCAT_VERSION
from pgxrecord.caller import GeneCall

_SCHEMA = """
CREATE TABLE IF NOT EXISTS records (
    record_id INTEGER PRIMARY KEY AUTOINCREMENT,
    subject_id TEXT NOT NULL,
    pharmcat_version TEXT NOT NULL,
    guideline_version TEXT NOT NULL,
    ingested_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS gene_calls (
    record_id INTEGER NOT NULL REFERENCES records(record_id),
    gene TEXT NOT NULL,
    diplotype TEXT,
    phenotype TEXT,
    coverage TEXT NOT NULL CHECK (
        coverage IN ('called', 'not_covered', 'indeterminate')
    ),
    PRIMARY KEY (record_id, gene)
);

CREATE INDEX IF NOT EXISTS idx_gene_calls_gene ON gene_calls(gene);
CREATE INDEX IF NOT EXISTS idx_records_subject ON records(subject_id);

CREATE TRIGGER IF NOT EXISTS gene_calls_no_update
BEFORE UPDATE ON gene_calls
BEGIN
    SELECT RAISE(ABORT, 'gene_calls is append-only');
END;

CREATE TRIGGER IF NOT EXISTS gene_calls_no_delete
BEFORE DELETE ON gene_calls
BEGIN
    SELECT RAISE(ABORT, 'gene_calls is append-only');
END;

CREATE TRIGGER IF NOT EXISTS records_no_update
BEFORE UPDATE ON records
BEGIN
    SELECT RAISE(ABORT, 'records is append-only');
END;

CREATE TRIGGER IF NOT EXISTS records_no_delete
BEFORE DELETE ON records
BEGIN
    SELECT RAISE(ABORT, 'records is append-only');
END;
"""


class RecordStore:
    """Append-only store of pharmacogenomic records."""

    def __init__(self, db_path: Path):
        self.db_path = db_path
        with self._connect() as conn:
            conn.executescript(_SCHEMA)

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        """Open a connection, commit on success, and always close.

        sqlite3's own connection context manager commits the transaction but
        does NOT close the handle, so using it directly leaks a file
        descriptor per call. This wraps both.
        """
        conn = sqlite3.connect(self.db_path)
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            with conn:
                yield conn
        finally:
            conn.close()

    def append(
        self,
        subject_id: str,
        calls: list[GeneCall],
        guideline_version: str,
    ) -> int:
        """Write a new record. Never modifies an existing one."""
        with self._connect() as conn:
            cursor = conn.execute(
                "INSERT INTO records "
                "(subject_id, pharmcat_version, guideline_version, ingested_at) "
                "VALUES (?, ?, ?, ?)",
                (
                    subject_id,
                    PHARMCAT_VERSION,
                    guideline_version,
                    datetime.now(timezone.utc).isoformat(),
                ),
            )
            record_id = int(cursor.lastrowid)
            conn.executemany(
                "INSERT INTO gene_calls "
                "(record_id, gene, diplotype, phenotype, coverage) "
                "VALUES (?, ?, ?, ?, ?)",
                [
                    (record_id, c.gene, c.diplotype, c.phenotype, c.coverage)
                    for c in calls
                ],
            )
        return record_id

    def history(self, subject_id: str) -> list[int]:
        """Every record_id for a subject, oldest first."""
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT record_id FROM records WHERE subject_id = ? "
                "ORDER BY record_id",
                (subject_id,),
            ).fetchall()
        return [row[0] for row in rows]

    def latest(self, subject_id: str) -> list[GeneCall]:
        """Gene calls from the subject's most recent record."""
        ids = self.history(subject_id)
        if not ids:
            return []
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT gene, diplotype, phenotype, coverage FROM gene_calls "
                "WHERE record_id = ? ORDER BY gene",
                (ids[-1],),
            ).fetchall()
        return [GeneCall(*row) for row in rows]

    def record_versions(self, record_id: int) -> tuple[str, str]:
        """The (pharmcat_version, guideline_version) a record was made with."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT pharmcat_version, guideline_version FROM records "
                "WHERE record_id = ?",
                (record_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"no record {record_id}")
        return (row[0], row[1])

    def subjects_with_gene(self, gene: str) -> list[str]:
        """Subjects with any stored call for this gene, across all records.

        Deliberately not restricted to the latest record: a guideline change
        is relevant to anyone ever genotyped for the gene.
        """
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT DISTINCT r.subject_id FROM records r "
                "JOIN gene_calls g ON g.record_id = r.record_id "
                "WHERE g.gene = ? ORDER BY r.subject_id",
                (gene,),
            ).fetchall()
        return [row[0] for row in rows]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_store.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pgxrecord/store.py tests/test_store.py
git commit -m "feat: append-only record store with schema-enforced immutability"
```

---

### Task 7: Drug query with the three-outcome invariant

**Files:**
- Create: `src/pgxrecord/guidelines.py`
- Create: `data/gene_drug_pairs.json`
- Create: `src/pgxrecord/evaluate.py`
- Test: `tests/test_evaluate.py`

**Interfaces:**
- Consumes: `GeneCall` (Task 5), `RecordStore` (Task 6)
- Produces:
  - `class GuidelineRef` — frozen dataclass, fields `gene: str`, `drug: str`, `cpic_pair_id: str`, `url: str`
  - `def load_pairs(path: Path) -> list[GuidelineRef]`
  - `def find_pairs_for_drug(drug: str, pairs: list[GuidelineRef]) -> list[GuidelineRef]`
  - `class QueryResult` — frozen dataclass, fields `outcome: str`, `gene: str | None`, `phenotype: str | None`, `guideline: GuidelineRef | None`, `explanation: str`
  - `def query_drug(store: RecordStore, subject_id: str, drug: str, pairs: list[GuidelineRef]) -> list[QueryResult]`

`outcome` is exactly one of `"guidance_found"`, `"no_guidance_for_pair"`, `"cannot_assess"`.

- [ ] **Step 1: Create the gene-drug pair reference**

`data/gene_drug_pairs.json` — identifiers and URLs only, no guideline prose:

```json
[
  {
    "gene": "CYP2C19",
    "drug": "clopidogrel",
    "cpic_pair_id": "CYP2C19-clopidogrel",
    "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
  },
  {
    "gene": "CYP2D6",
    "drug": "codeine",
    "cpic_pair_id": "CYP2D6-codeine",
    "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
  },
  {
    "gene": "DPYD",
    "drug": "fluorouracil",
    "cpic_pair_id": "DPYD-fluorouracil",
    "url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/"
  },
  {
    "gene": "TPMT",
    "drug": "azathioprine",
    "cpic_pair_id": "TPMT-azathioprine",
    "url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt-and-nudt15/"
  },
  {
    "gene": "SLCO1B1",
    "drug": "simvastatin",
    "cpic_pair_id": "SLCO1B1-simvastatin",
    "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/"
  }
]
```

- [ ] **Step 2: Write the failing test**

`tests/test_evaluate.py`:

```python
from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall
from pgxrecord.evaluate import query_drug
from pgxrecord.guidelines import load_pairs
from pgxrecord.store import RecordStore

PAIRS = load_pairs(Path(__file__).resolve().parents[1] / "data/gene_drug_pairs.json")


@pytest.fixture
def store(tmp_path):
    return RecordStore(tmp_path / "records.db")


def test_guidance_found_for_called_gene(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "clopidogrel", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "guidance_found"
    assert results[0].phenotype == "Intermediate Metabolizer"
    assert results[0].guideline.cpic_pair_id == "CYP2C19-clopidogrel"
    assert results[0].guideline.url.startswith("https://")


def test_no_guidance_for_drug_with_no_cpic_pair(store):
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "amoxicillin", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "no_guidance_for_pair"
    assert results[0].guideline is None


def test_cannot_assess_when_gene_not_covered(store):
    """The invariant. A gene the array never informed carries NO information."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)

    assert len(results) == 1
    assert results[0].outcome == "cannot_assess"
    assert results[0].phenotype is None


def test_cannot_assess_when_gene_absent_from_record(store):
    """A gene missing entirely is unassessable, never 'no interaction'."""
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*1", "Normal Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)

    assert [r.outcome for r in results] == ["cannot_assess"]


def test_cannot_assess_when_indeterminate(store):
    store.append(
        "s1",
        [GeneCall("TPMT", None, None, "indeterminate")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "azathioprine", PAIRS)

    assert results[0].outcome == "cannot_assess"


def test_cannot_assess_is_never_worded_as_reassurance(store):
    """No 'cannot_assess' explanation may read as an all-clear."""
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "codeine", PAIRS)
    text = results[0].explanation.lower()

    for reassuring in (
        "no interaction",
        "no issue",
        "safe",
        "normal",
        "no guidance",
        "clear",
        "fine",
    ):
        assert reassuring not in text, f"reassuring phrase {reassuring!r} in {text!r}"
    assert "cannot" in text or "unable" in text


def test_no_output_contains_a_dose_or_clinical_directive(store):
    """No clinical claims. Guidance is referenced, never restated as advice."""
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    results = query_drug(store, "s1", "clopidogrel", PAIRS)
    text = results[0].explanation.lower()

    for directive in ("mg", "you should", "take ", "avoid ", "dose of", "recommend"):
        assert directive not in text, f"clinical directive {directive!r} in {text!r}"


def test_unknown_subject_cannot_be_assessed(store):
    results = query_drug(store, "nobody", "clopidogrel", PAIRS)
    assert results[0].outcome == "cannot_assess"


def test_every_outcome_is_one_of_three_values(store):
    store.append(
        "s1",
        [
            GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called"),
            GeneCall("CYP2D6", None, None, "not_covered"),
        ],
        guideline_version="cpic-2026-07",
    )

    outcomes = {
        r.outcome
        for drug in ("clopidogrel", "codeine", "amoxicillin")
        for r in query_drug(store, "s1", drug, PAIRS)
    }
    assert outcomes <= {"guidance_found", "no_guidance_for_pair", "cannot_assess"}
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_evaluate.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.evaluate'`

- [ ] **Step 4: Write the guideline reference loader**

`src/pgxrecord/guidelines.py`:

```python
"""Gene-drug pair references.

We store identifiers and URLs only, never guideline prose. Two reasons:
guidelines are revised, so a stored copy is a stale copy; and redistributing
CPIC/PharmGKB content commercially may require permission we have not
obtained. Linking sidesteps both.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GuidelineRef:
    """A pointer to published guidance for one gene-drug pair."""

    gene: str
    drug: str
    cpic_pair_id: str
    url: str


def load_pairs(path: Path) -> list[GuidelineRef]:
    """Load the gene-drug pair reference table."""
    return [GuidelineRef(**entry) for entry in json.loads(path.read_text())]


def find_pairs_for_drug(drug: str, pairs: list[GuidelineRef]) -> list[GuidelineRef]:
    """Every gene-drug pair matching a drug name, case-insensitively."""
    needle = drug.strip().lower()
    return [p for p in pairs if p.drug.lower() == needle]
```

- [ ] **Step 5: Write the evaluator**

`src/pgxrecord/evaluate.py`:

```python
"""Answer drug queries against a stored record.

The load-bearing rule of this module: absence of guidance and absence of data
are different answers and must never collapse into one. A gene a consumer
array never covered carries NO information, and presenting that as "no
interaction found" is the most dangerous thing this system could do.

Hence three outcomes, always distinguishable:

  guidance_found       gene called, CPIC publishes guidance for this pair
  no_guidance_for_pair gene called, CPIC publishes nothing for this drug
  cannot_assess        gene not covered, indeterminate, or absent -- unknown
"""

from __future__ import annotations

from dataclasses import dataclass

from pgxrecord.caller import CALLED
from pgxrecord.guidelines import GuidelineRef, find_pairs_for_drug
from pgxrecord.store import RecordStore

GUIDANCE_FOUND = "guidance_found"
NO_GUIDANCE_FOR_PAIR = "no_guidance_for_pair"
CANNOT_ASSESS = "cannot_assess"


@dataclass(frozen=True)
class QueryResult:
    """One gene's answer for a drug query."""

    outcome: str
    gene: str | None
    phenotype: str | None
    guideline: GuidelineRef | None
    explanation: str


def query_drug(
    store: RecordStore,
    subject_id: str,
    drug: str,
    pairs: list[GuidelineRef],
) -> list[QueryResult]:
    """Return one result per gene-drug pair relevant to this drug."""
    relevant = find_pairs_for_drug(drug, pairs)
    if not relevant:
        return [
            QueryResult(
                outcome=NO_GUIDANCE_FOR_PAIR,
                gene=None,
                phenotype=None,
                guideline=None,
                explanation=(
                    f"CPIC publishes no gene-drug pair for {drug!r} in this "
                    f"reference table."
                ),
            )
        ]

    calls = {c.gene: c for c in store.latest(subject_id)}
    results: list[QueryResult] = []

    for pair in relevant:
        call = calls.get(pair.gene)

        if call is None:
            results.append(
                QueryResult(
                    outcome=CANNOT_ASSESS,
                    gene=pair.gene,
                    phenotype=None,
                    guideline=pair,
                    explanation=(
                        f"Cannot assess {pair.gene}: this subject has no stored "
                        f"call for it, so the genotype is unknown."
                    ),
                )
            )
            continue

        if call.coverage != CALLED:
            results.append(
                QueryResult(
                    outcome=CANNOT_ASSESS,
                    gene=pair.gene,
                    phenotype=None,
                    guideline=pair,
                    explanation=(
                        f"Cannot assess {pair.gene}: coverage is "
                        f"{call.coverage!r}, so the genotype is unknown. This is "
                        f"absence of data, not absence of an interaction."
                    ),
                )
            )
            continue

        results.append(
            QueryResult(
                outcome=GUIDANCE_FOUND,
                gene=pair.gene,
                phenotype=call.phenotype,
                guideline=pair,
                explanation=(
                    f"{pair.gene} genotype {call.diplotype} corresponds to "
                    f"phenotype {call.phenotype!r}. CPIC publishes guidance for "
                    f"pair {pair.cpic_pair_id}; see {pair.url} for the "
                    f"guideline text."
                ),
            )
        )

    return results
```

- [ ] **Step 6: Run test to verify it passes**

Run: `pytest tests/test_evaluate.py -v`
Expected: PASS (9 tests)

- [ ] **Step 7: Commit**

```bash
git add src/pgxrecord/guidelines.py src/pgxrecord/evaluate.py data/gene_drug_pairs.json tests/test_evaluate.py
git commit -m "feat: drug query with three-outcome invariant and linked guidelines"
```

---

### Task 8: Guideline-bump diffing

**Files:**
- Create: `src/pgxrecord/drift.py`
- Test: `tests/test_drift.py`

**Interfaces:**
- Consumes: `RecordStore` (Task 6), `GuidelineRef` (Task 7)
- Produces:
  - `class AffectedRecord` — frozen dataclass, fields `subject_id: str`, `gene: str`, `changed_pair_ids: list[str]`
  - `def affected_by_guideline_change(store: RecordStore, changed_pair_ids: set[str], pairs: list[GuidelineRef]) -> list[AffectedRecord]`

This is the compounding behavior a stateless tool cannot offer: stored records gain value when guidelines are revised.

- [ ] **Step 1: Write the failing test**

`tests/test_drift.py`:

```python
from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall
from pgxrecord.drift import affected_by_guideline_change
from pgxrecord.guidelines import load_pairs
from pgxrecord.store import RecordStore

PAIRS = load_pairs(Path(__file__).resolve().parents[1] / "data/gene_drug_pairs.json")


@pytest.fixture
def store(tmp_path):
    store = RecordStore(tmp_path / "records.db")
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )
    store.append(
        "s2",
        [GeneCall("DPYD", "Ref/Ref", "Normal Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )
    store.append(
        "s3",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )
    return store


def test_finds_subjects_affected_by_a_changed_pair(store):
    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel"}, PAIRS
    )

    assert len(affected) == 1
    assert affected[0].subject_id == "s1"
    assert affected[0].gene == "CYP2C19"
    assert affected[0].changed_pair_ids == ["CYP2C19-clopidogrel"]


def test_unaffected_subjects_are_not_reported(store):
    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel"}, PAIRS
    )
    assert {a.subject_id for a in affected} == {"s1"}


def test_multiple_changed_pairs(store):
    affected = affected_by_guideline_change(
        store, {"CYP2C19-clopidogrel", "DPYD-fluorouracil"}, PAIRS
    )
    assert {a.subject_id for a in affected} == {"s1", "s2"}


def test_not_covered_genes_are_still_reported(store):
    """A not_covered gene still matters: new guidance may warrant re-testing."""
    affected = affected_by_guideline_change(store, {"CYP2D6-codeine"}, PAIRS)

    assert len(affected) == 1
    assert affected[0].subject_id == "s3"
    assert affected[0].gene == "CYP2D6"


def test_unknown_pair_id_affects_nobody(store):
    assert affected_by_guideline_change(store, {"NOPE-nothing"}, PAIRS) == []


def test_no_changes_affects_nobody(store):
    assert affected_by_guideline_change(store, set(), PAIRS) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_drift.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.drift'`

- [ ] **Step 3: Write the implementation**

`src/pgxrecord/drift.py`:

```python
"""Report which stored records a guideline revision touches.

This is the reason the record store exists. A batch tool answers "what does
this genotype mean today". A persistent record can answer "whose stored
results changed meaning when the guidance moved" -- including subjects whose
gene was never covered, since new guidance may justify proper testing.
"""

from __future__ import annotations

from dataclasses import dataclass

from pgxrecord.guidelines import GuidelineRef
from pgxrecord.store import RecordStore


@dataclass(frozen=True)
class AffectedRecord:
    """A stored subject-gene pair touched by a guideline change."""

    subject_id: str
    gene: str
    changed_pair_ids: list[str]


def affected_by_guideline_change(
    store: RecordStore,
    changed_pair_ids: set[str],
    pairs: list[GuidelineRef],
) -> list[AffectedRecord]:
    """Find stored records involving gene-drug pairs whose guidance changed."""
    if not changed_pair_ids:
        return []

    genes_to_pairs: dict[str, list[str]] = {}
    for pair in pairs:
        if pair.cpic_pair_id in changed_pair_ids:
            genes_to_pairs.setdefault(pair.gene, []).append(pair.cpic_pair_id)

    affected: list[AffectedRecord] = []
    for gene, pair_ids in sorted(genes_to_pairs.items()):
        for subject_id in store.subjects_with_gene(gene):
            affected.append(
                AffectedRecord(
                    subject_id=subject_id,
                    gene=gene,
                    changed_pair_ids=sorted(pair_ids),
                )
            )
    return affected
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_drift.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pgxrecord/drift.py tests/test_drift.py
git commit -m "feat: report stored records affected by guideline revisions"
```

---

### Task 9: End-to-end CLI and README

**Files:**
- Create: `src/pgxrecord/cli.py`
- Modify: `pyproject.toml` (add `[project.scripts]`)
- Create: `README.md`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: every prior task
- Produces: console script `pgxrecord` with subcommands `ingest` and `query`

- [ ] **Step 1: Write the failing test**

`tests/test_cli.py`:

```python
import json
from pathlib import Path

import pytest

from pgxrecord.caller import GeneCall
from pgxrecord.cli import cmd_query, ingest_to_calls
from pgxrecord.store import RecordStore

FIXTURES = Path("tests/fixtures")


def test_ingest_to_calls_reports_uncovered_genes_without_docker(tmp_path):
    """The ingest half runs without Docker and must surface coverage."""
    vcf_path, report = ingest_to_calls(
        FIXTURES / "23andme_valid.txt",
        Path(__file__).resolve().parents[1] / "data/pharmcat_positions_3.4.0.vcf",
        tmp_path,
    )

    assert vcf_path.is_file()
    # The fixture covers only a couple of DPYD/CYP2C19 positions, so CYP2D6
    # -- which needs rsID-less positions -- must be fully uncovered.
    assert "CYP2D6" in report.genes_fully_uncovered


def test_query_prints_cannot_assess_prominently(tmp_path, capsys):
    store = RecordStore(tmp_path / "records.db")
    store.append(
        "s1",
        [GeneCall("CYP2D6", None, None, "not_covered")],
        guideline_version="cpic-2026-07",
    )

    cmd_query(store, "s1", "codeine")
    out = capsys.readouterr().out

    assert "CANNOT ASSESS" in out
    assert "no interaction" not in out.lower()


def test_query_prints_guidance_with_citation(tmp_path, capsys):
    store = RecordStore(tmp_path / "records.db")
    store.append(
        "s1",
        [GeneCall("CYP2C19", "*1/*2", "Intermediate Metabolizer", "called")],
        guideline_version="cpic-2026-07",
    )

    cmd_query(store, "s1", "clopidogrel")
    out = capsys.readouterr().out

    assert "CYP2C19-clopidogrel" in out
    assert "https://" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cli.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgxrecord.cli'`

- [ ] **Step 3: Write the CLI**

`src/pgxrecord/cli.py`:

```python
"""Command line interface.

Output wording is constrained: cannot_assess is printed as a loud CANNOT
ASSESS banner, never as a quiet empty result. See evaluate.py for why.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pgxrecord import POSITIONS_FILENAME
from pgxrecord.caller import parse_phenotype_json, run_pharmcat
from pgxrecord.evaluate import CANNOT_ASSESS, GUIDANCE_FOUND, query_drug
from pgxrecord.guidelines import load_pairs
from pgxrecord.ingest.raw import parse_23andme
from pgxrecord.ingest.vcf import CoverageReport, build_vcf
from pgxrecord.positions import load_positions
from pgxrecord.store import RecordStore

# Anchored to the installed package, not the CWD, so the CLI works from
# anywhere rather than only from the repo root.
DATA_DIR = Path(__file__).resolve().parents[2].parent / "data"
PAIRS_PATH = DATA_DIR / "gene_drug_pairs.json"


def ingest_to_calls(
    raw_path: Path, positions_path: Path, workdir: Path
) -> tuple[Path, CoverageReport]:
    """Parse a raw file and write a VCF. Does not require Docker."""
    calls = parse_23andme(raw_path)
    positions = load_positions(positions_path)
    vcf_path = workdir / f"{raw_path.stem}.vcf"
    report = build_vcf(calls, positions, vcf_path)
    return vcf_path, report


def cmd_ingest(
    store: RecordStore, raw_path: Path, subject_id: str, workdir: Path
) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    vcf_path, report = ingest_to_calls(
        raw_path, DATA_DIR / POSITIONS_FILENAME, workdir
    )
    print(f"wrote {vcf_path}")
    print(f"covered positions: {len(report.covered_rsids)}")
    print(f"uncovered positions: {len(report.uncovered_rsids)}")
    print(f"genes with no coverage: {', '.join(sorted(report.genes_fully_uncovered))}")

    phenotype_json = run_pharmcat(vcf_path, workdir)
    gene_calls = parse_phenotype_json(
        phenotype_json, uncovered_genes=report.genes_fully_uncovered
    )
    record_id = store.append(subject_id, gene_calls, guideline_version="cpic-2026-07")
    print(f"stored record {record_id} for subject {subject_id}")


def cmd_query(store: RecordStore, subject_id: str, drug: str) -> None:
    pairs = load_pairs(PAIRS_PATH)
    for result in query_drug(store, subject_id, drug, pairs):
        if result.outcome == CANNOT_ASSESS:
            print(f"[CANNOT ASSESS] {result.explanation}")
        elif result.outcome == GUIDANCE_FOUND:
            print(f"[GUIDANCE] {result.explanation}")
        else:
            print(f"[NO CPIC PAIR] {result.explanation}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pgxrecord")
    parser.add_argument("--db", type=Path, default=Path("records.db"))
    sub = parser.add_subparsers(dest="command", required=True)

    p_ingest = sub.add_parser("ingest", help="ingest a 23andMe raw file")
    p_ingest.add_argument("raw_path", type=Path)
    p_ingest.add_argument("--subject", required=True)
    p_ingest.add_argument("--workdir", type=Path, default=Path("work"))

    p_query = sub.add_parser("query", help="query a drug against a stored record")
    p_query.add_argument("drug")
    p_query.add_argument("--subject", required=True)

    args = parser.parse_args(argv)
    store = RecordStore(args.db)

    if args.command == "ingest":
        cmd_ingest(store, args.raw_path, args.subject, args.workdir)
    else:
        cmd_query(store, args.subject, args.drug)
    return 0
```

Add to `pyproject.toml`:

```toml
[project.scripts]
pgxrecord = "pgxrecord.cli:main"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pip install -e ".[dev]" && pytest tests/test_cli.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Write the README**

`README.md`:

```markdown
# pgxrecord

A longitudinal pharmacogenomic record built on [PharmCAT](https://github.com/PharmGKB/PharmCAT).

## What this is

PharmCAT is excellent and stateless: VCF in, report out. This project adds the
layer it does not have — a persistent, version-stamped record that is written
once and re-evaluated every time a new drug is queried or a guideline is revised.

**This is research and reference tooling. It is not clinical decision support,
it makes no claim about any person's care, and it must not be used to make
treatment decisions.**

## What it does not do

- Call star alleles. PharmCAT does that; we pin `pgkb/pharmcat:3.4.0` and never reimplement it.
- Restate guideline text. We store CPIC gene-drug pair identifiers and link out.
- Resolve `CYP2D6`. That gene depends on copy-number and structural variation
  that consumer arrays cannot resolve. It is reported `not_covered`.

## Coverage honesty

Consumer arrays genotype a sparse subset of positions. The system distinguishes
three states and never collapses them:

| State | Meaning |
|---|---|
| `guidance_found` | Gene called, CPIC publishes guidance for this pair |
| `no_guidance_for_pair` | Gene called, CPIC publishes nothing for this drug |
| `cannot_assess` | Gene not covered or indeterminate — **we do not know** |

`cannot_assess` is never rendered as reassurance. Absence of data is not
absence of an interaction.

## Usage

```bash
./scripts/fetch_positions.sh
pip install -e ".[dev]"

pgxrecord ingest genome.txt --subject me     # requires Docker
pgxrecord query clopidogrel --subject me
```

## Limitations

- Consumer genotype data is not clinically confirmed.
- Guidance is pinned to a version; the tool is only as current as that pin.
- Star-allele frequency and guideline evidence are unevenly distributed across
  ancestries, so guidance quality is not uniform.

## License and attribution

PharmCAT is MPL-2.0. CPIC and PharmGKB guideline content is referenced by link,
not redistributed. Embedding guideline text or commercializing this project
requires resolving CPIC/PharmGKB data-use terms first.
```

- [ ] **Step 6: Run the full suite**

Run: `pytest -v`
Expected: PASS, 53 tests across 8 files

- [ ] **Step 7: Commit**

```bash
git add src/pgxrecord/cli.py pyproject.toml README.md tests/test_cli.py
git commit -m "feat: CLI for ingest and drug query, plus README"
```

---

## Self-Review

**Spec coverage**

| Spec section | Task |
|---|---|
| Ingest (vendor/build detection, rsID join) | 3, 4 |
| Call (pinned PharmCAT wrapper) | 5 |
| Record (append-only, version-stamped) | 6 |
| Re-evaluate (drug query) | 7 |
| Re-evaluate (guideline bump) | 8 |
| Three coverage states | 4, 5 |
| Three query outcomes + invariant | 7 |
| Error handling (reject, never guess) | 3, 5 |
| Testing obligations | every task |
| Known limitations documented | 9 (README) |
| Guidelines linked not embedded | 7 |

No spec requirement is unimplemented. The spec's "clinical VCF ingest" and "FHIR integration" are explicitly out of v1 scope and correctly absent.

**Placeholder scan:** none. Every step contains runnable code or an exact command.

**Type consistency:** `GeneCall` fields (`gene`, `diplotype`, `phenotype`, `coverage`) are identical across Tasks 5–9. `CoverageReport.genes_fully_uncovered` feeds `parse_phenotype_json(uncovered_genes=...)` in Tasks 5 and 9 under the same name. Coverage strings `called`/`not_covered`/`indeterminate` match the SQLite `CHECK` constraint in Task 6 exactly. `GuidelineRef` field names match `data/gene_drug_pairs.json` keys, which is required because `load_pairs` uses `GuidelineRef(**entry)`.

**Note on Task 9:** `test_ingest_to_calls_reports_uncovered_genes_without_docker` exercises ingest only. Full end-to-end ingest requires Docker and is not covered by an automated test — run it manually with a real 23andMe export.
