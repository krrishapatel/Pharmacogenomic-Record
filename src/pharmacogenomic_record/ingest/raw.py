"""Parse 23andMe raw genotype exports.

Consumer array files are messy in specific known ways, and this module is
where that mess is contained. Two rules:

1. Never guess the vendor or genome build. An unrecognized header is a
   rejection, not a default.
2. Drop rows that cannot be joined -- internal 'i' identifiers and no-calls
   ('--'). Dropping them here means downstream code sees only usable calls,
   and the positions they would have covered surface as not_covered.

Hemizygous single-allele calls are the exception to rule 2: they are kept and
passed on, even though nothing downstream can encode them yet. They used to be
dropped here on the stated grounds that no valid VCF genotype existed for them,
which is not true -- a haploid GT of "0" or "1" is valid VCF 4.2 and is the
standard encoding for male non-PAR chrX, chrY and chrM. Dropping them silently
made a measured position indistinguishable from one the array never had, and for
a male sample that is every G6PD position in the reference table: the only
X-linked gene there, and the one where hemizygosity is the clinical point. They
are carried through so the coverage report can say so out loud.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# Genotypes that carry no information. A single dash is a hemizygous no-call
# (chrY/MT); a double dash is the diploid form.
NO_CALLS = frozenset({"--", "-"})

# 23andMe declares the build in a fixed sentence. We capture the number and
# require it to be 37 rather than searching for the substring "build 37"
# anywhere in the header -- a build-38 file whose header merely mentions
# build 37 in prose would otherwise pass, and every coordinate downstream
# would be silently wrong.
_BUILD_RE = re.compile(r"reference human assembly build (\d+)")

# Indel calls (D=deletion, I=insertion) are not nucleotide alleles and cannot
# be matched against the reference ref/alt bases, so they are unjoinable.
_INDEL_CODES = frozenset({"D", "I"})


class UnsupportedRawFile(Exception):
    """The raw file is not a format we can convert safely."""


@dataclass(frozen=True)
class RawCall:
    """One genotype call from a consumer array."""

    rsid: str
    chrom: str
    pos: int
    genotype: str

    @property
    def is_hemizygous(self) -> bool:
        """True for a single-allele call, as reported for male non-PAR chrX.

        A property rather than a stored flag so it cannot disagree with the
        genotype it describes. No-calls and indel codes never reach a RawCall,
        so a one-character genotype here is a real single nucleotide.
        """
        return len(self.genotype) == 1


def _validate_header(header_text: str) -> None:
    lowered = header_text.lower()
    if "23andme" not in lowered:
        raise UnsupportedRawFile(
            "no recognizable 23andMe header found; refusing to guess the "
            "vendor or genome build"
        )

    builds = set(_BUILD_RE.findall(lowered))
    if not builds:
        raise UnsupportedRawFile(
            "header declares no reference assembly build; refusing to guess "
            "the genome build"
        )
    if len(builds) > 1:
        # Taking the first match here would be a coin flip on whether every
        # coordinate downstream is right.
        raise UnsupportedRawFile(
            f"header declares more than one reference assembly build "
            f"({', '.join(sorted(builds))}); refusing to guess which applies"
        )
    build = builds.pop()
    if build != "37":
        raise UnsupportedRawFile(
            f"expected reference assembly build 37, header declares build "
            f"{build}; this would invalidate the rsID join"
        )


def parse_23andme(path: Path) -> list[RawCall]:
    """Parse a 23andMe raw export into joinable genotype calls.

    Rows that carry no usable information are dropped: internal 'i'
    identifiers, no-calls, and indel codes. Rows whose *shape* is wrong are a
    rejection, not a skip -- a wrong-format file must not parse to an empty list
    and look like a clean file with nothing relevant in it.

    Hemizygous single-allele calls are kept, flagged by `RawCall.is_hemizygous`.
    Nothing writes them to a VCF yet, but a measured position has to stay
    distinguishable from an absent one; see the module docstring.
    """
    lines = path.read_text().splitlines()
    header_text = "\n".join(line for line in lines if line.startswith("#"))
    _validate_header(header_text)

    calls: list[RawCall] = []
    for number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) != 4:
            raise UnsupportedRawFile(
                f"{path}:{number}: expected 4 tab-separated columns, found "
                f"{len(fields)}; this does not look like a 23andMe export"
            )

        rsid, chrom, pos, genotype = fields
        if not pos.isdigit():
            raise UnsupportedRawFile(
                f"{path}:{number}: position {pos!r} is not a number"
            )

        if not rsid.startswith("rs"):
            continue
        if genotype in NO_CALLS:
            continue
        if set(genotype) & _INDEL_CODES:
            continue
        if len(genotype) not in (1, 2):
            # Neither haploid nor diploid, so it is not a genotype this format
            # is supposed to contain.
            continue

        calls.append(
            RawCall(rsid=rsid, chrom=chrom, pos=int(pos), genotype=genotype)
        )

    if not calls:
        raise UnsupportedRawFile(
            f"{path}: no usable genotype calls found; the file parsed but "
            f"produced nothing, which usually means the format is wrong"
        )
    return calls
