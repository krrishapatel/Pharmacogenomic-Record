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
