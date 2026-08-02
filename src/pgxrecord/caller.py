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
from collections.abc import Iterable
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


def parse_phenotype_json(
    path: Path, uncovered_genes: Iterable[str]
) -> list[GeneCall]:
    """Translate PharmCAT phenotype output into GeneCall records.

    uncovered_genes always wins. If the array never informed a gene, no
    PharmCAT output about it can be treated as a call. Accepts any iterable of
    gene names (Task 9 passes a frozenset from CoverageReport).
    """
    uncovered = frozenset(uncovered_genes)

    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as err:
        raise PharmcatError(f"could not parse PharmCAT output {path}: {err}") from err

    calls: list[GeneCall] = []
    seen: set[str] = set()

    for gene, entry in (payload.get("phenotypes") or {}).items():
        seen.add(gene)
        if gene in uncovered:
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

    for gene in sorted(uncovered - seen):
        calls.append(GeneCall(gene, None, None, NOT_COVERED))

    return calls
