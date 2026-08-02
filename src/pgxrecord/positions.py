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
