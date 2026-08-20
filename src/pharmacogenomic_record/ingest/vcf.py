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

from pharmacogenomic_record.ingest.raw import RawCall
from pharmacogenomic_record.positions import ReferencePosition, index_by_rsid

_VCF_META = """##fileformat=VCFv4.2
##source=pharmacogenomic_record
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

    The three gene sets partition every gene in the reference table, and the
    boundary between them is the STRICT rule: a gene is eligible to be called
    only when EVERY rsID-joinable position of it was covered.

      genes_fully_covered      every joinable position covered -> may be called
      genes_partially_covered  some but not all covered -> indeterminate
      genes_fully_uncovered    no position covered at all -> not_covered

    Strict, with no threshold and no ratio, because there is no principled one:
    PharmCAT assumes reference at any position it was not given, so a gene with
    39 of 40 positions covered still yields a confident "*1/*1 Normal
    Metabolizer", and the one missing position is exactly where a variant would
    have been. Any cutoff below "all of them" is a number nobody can defend.

    `genes_fully_uncovered` is decided FIRST, and that ordering is load-bearing.
    A gene whose every position lacks an rsID has an empty joinable set, so
    "every joinable position covered" is vacuously true for it -- ranked the
    other way round, a gene the array said nothing whatsoever about would come
    back eligible to be called. Membership here requires at least one actually
    covered position, so the vacuous case lands in fully_uncovered where it
    belongs.

    What this rule does NOT close: rsID-less positions are excluded from the
    denominator, so a gene can be fully covered while real positions of it were
    never observed -- G6PD is only 39% rsID-bearing, RYR1 78%, and 9 of the 22
    genes hold at least one such position. `unjoinable_positions` is reported
    for exactly that reason; it is a known residual gap, not a solved one.

    Fields are frozensets: `frozen=True` stops rebinding but not
    `report.covered_rsids.add(...)`, and a mutable set field would also make
    the report unhashable and unsafe to share.

    `unjoinable_positions` is the count of reference positions that have no
    rsID at all. They appear in neither covered nor uncovered because they have
    no rsID to key on -- but they are real gaps in coverage, so the count is
    surfaced rather than silently dropped.

    `hemizygous_rsids` and `hemizygous_genes` are the same idea for the same
    reason: reference positions the array did measure, as a single allele, that
    this pipeline cannot yet write to a VCF. They count as uncovered, because
    nothing about them reaches PharmCAT, but "measured and not encodable" is a
    different fact from "never measured" and the two must not be reported as one.
    For a male sample this is every joinable G6PD position, G6PD being the only
    X-linked gene in the table.
    """

    covered_rsids: frozenset[str]
    uncovered_rsids: frozenset[str]
    genes_fully_uncovered: frozenset[str]
    genes_partially_covered: frozenset[str]
    genes_fully_covered: frozenset[str]
    unjoinable_positions: int
    hemizygous_rsids: frozenset[str]
    hemizygous_genes: frozenset[str]


def translate_genotype(genotype: str, ref: ReferencePosition) -> str | None:
    """Convert raw allele letters to a VCF numeric genotype.

    Returns None when any allele is neither the reference nor a known
    alternate, which means the call tells us nothing about this position.

    A single-allele hemizygous call also returns None. That is a limitation, not
    a property of VCF: a haploid GT is legal and is what such a call should
    become. It is left unencoded because whether the pinned PharmCAT image
    accepts a haploid GT for G6PD has not been verified here, and sending it
    something unverified is worse than reporting the gap. `CoverageReport`
    reports the gap.
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
    genes_with_any_coverage = {
        by_rsid[r].gene for r in covered if by_rsid[r].gene is not None
    }
    all_genes = {p.gene for p in positions if p.gene is not None}

    # The denominator of the strict rule, per gene: every rsID-joinable position
    # the reference table lists for it. Built from `positions` rather than from
    # `by_rsid.values()` so that a future duplicate rsID -- which the dict would
    # silently collapse to one entry -- cannot shrink a gene's denominator and
    # promote a partially covered gene to fully covered. There are no duplicates
    # in 3.4.0; the point is that this does not depend on that staying true.
    joinable_by_gene: dict[str, set[str]] = {}
    for position in positions:
        if position.gene is not None and position.rsid is not None:
            joinable_by_gene.setdefault(position.gene, set()).add(position.rsid)

    # Ranked, not independent: fully_uncovered first, so that a gene with no
    # joinable positions at all cannot satisfy "every joinable position covered"
    # vacuously. See CoverageReport.
    genes_fully_uncovered = all_genes - genes_with_any_coverage
    genes_fully_covered = {
        gene
        for gene in genes_with_any_coverage
        if joinable_by_gene.get(gene, set()) <= covered
    }

    # Measured as a single allele at a position the reference table needs. Not
    # covered, because no VCF row was written for it, but distinct from a
    # position the array never reported at all.
    hemizygous_rsids = {
        call.rsid
        for call in calls
        if call.is_hemizygous and call.rsid in by_rsid
    } - covered
    hemizygous_genes = {
        by_rsid[rsid].gene
        for rsid in hemizygous_rsids
        if by_rsid[rsid].gene is not None
    }

    return CoverageReport(
        covered_rsids=frozenset(covered),
        uncovered_rsids=frozenset(all_joinable - covered),
        genes_fully_uncovered=frozenset(genes_fully_uncovered),
        genes_partially_covered=frozenset(
            genes_with_any_coverage - genes_fully_covered
        ),
        genes_fully_covered=frozenset(genes_fully_covered),
        unjoinable_positions=len(positions) - len(by_rsid),
        hemizygous_rsids=frozenset(hemizygous_rsids),
        hemizygous_genes=frozenset(hemizygous_genes),
    )
