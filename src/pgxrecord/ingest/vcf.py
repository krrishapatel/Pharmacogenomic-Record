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
    """Which reference positions the array actually informed.

    Fields are frozensets: `frozen=True` stops rebinding but not
    `report.covered_rsids.add(...)`, and a mutable set field would also make
    the report unhashable and unsafe to share.

    `unjoinable_positions` is the count of reference positions that have no
    rsID at all (208 in 3.4.0). They appear in neither covered nor uncovered
    because they have no rsID to key on -- but they are real gaps in
    coverage, so the count is surfaced rather than silently dropped.
    """

    covered_rsids: frozenset[str]
    uncovered_rsids: frozenset[str]
    genes_fully_uncovered: frozenset[str]
    genes_partially_covered: frozenset[str]
    unjoinable_positions: int


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
        covered_rsids=frozenset(covered),
        uncovered_rsids=frozenset(all_joinable - covered),
        genes_fully_uncovered=frozenset(all_genes - genes_covered_partly),
        genes_partially_covered=frozenset(genes_covered_partly),
        unjoinable_positions=len(positions) - len(by_rsid),
    )
