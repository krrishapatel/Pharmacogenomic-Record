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
import re
import shutil
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from pgxrecord import PHARMCAT_IMAGE

CALLED = "called"
NOT_COVERED = "not_covered"
INDETERMINATE = "indeterminate"

# Phenotype strings that mean "PharmCAT could not determine a phenotype".
# These are matched against the *phenotype* only, never against the diplotype
# label: labels have the form "allele1/allele2", so a substring marker
# containing "/" straddles the separator ("n/a" is a substring of
# "Canton/Aures", a real and clinically actionable G6PD-deficient result).
# Matching is whole-token, not raw substring, so "No Function" does not trip
# "no result".
_INDETERMINATE_PHENOTYPES = ("indeterminate", "unknown", "n/a", "no result")
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


def _vcf_basename(vcf_path: Path) -> str:
    """Strip .vcf / .vcf.gz to get PharmCAT's output base name.

    Path.stem removes exactly one suffix, so "sample.vcf.gz" would leave
    "sample.vcf" and we would look for "sample.vcf.phenotype.json".
    """
    name = vcf_path.name
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    if name.endswith(".vcf"):
        name = name[: -len(".vcf")]
    return name


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

    phenotype = workdir / f"{_vcf_basename(vcf_path)}.phenotype.json"
    if not phenotype.is_file():
        raise PharmcatError(
            f"PharmCAT produced no phenotype JSON at {phenotype}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return phenotype


def _tokens(text: str) -> list[str]:
    """Lowercase alphanumeric tokens, with every separator dropped.

    "N/A" -> ["n", "a"], "No Function" -> ["no", "function"].
    """
    return re.findall(r"[a-z0-9]+", text.lower())


_INDETERMINATE_MARKERS = tuple(
    _tokens(marker) for marker in _INDETERMINATE_PHENOTYPES
)


def _phenotype_is_indeterminate(phenotype: str | None) -> bool:
    """True if the phenotype string itself says "could not determine".

    An absent phenotype is NOT indeterminate: F2, F5, VKORC1, CFTR, IFNL3 and
    ABCG2 have no metabolizer phenotype at all, and absence of a metabolizer
    phenotype is not absence of a call. Markers are matched as contiguous
    whole-token runs, so "No Function" does not trip the "no result" marker.
    """
    if not phenotype:
        return False
    tokens = _tokens(phenotype)
    for marker in _INDETERMINATE_MARKERS:
        span = len(marker)
        if span and any(
            tokens[i:i + span] == marker for i in range(len(tokens) - span + 1)
        ):
            return True
    return False


def _allele_is_called(allele: object, gene: str) -> bool:
    """True if an allele slot holds a real, named allele.

    A null allele, or an allele whose name is null/missing/blank, is PharmCAT's
    representation of a no-call -- this is the structural signal, and it is the
    primary one. Anything that is neither null nor a dict is unparseable.
    """
    if allele is None:
        return False
    if not isinstance(allele, dict):
        raise PharmcatError(
            f"PharmCAT output for {gene} has a non-object allele: {allele!r}"
        )
    name = allele.get("name")
    if name is None:
        return False
    if not isinstance(name, str):
        raise PharmcatError(
            f"PharmCAT output for {gene} has a non-string allele name: {name!r}"
        )
    return bool(name.strip())


def _resolve(gene: str, entry: object) -> GeneCall:
    """Turn one PharmCAT phenotypes entry into a GeneCall.

    Classification order, strictest first:

    1. No candidate diplotypes at all -> indeterminate.
    2. More than one candidate with differing labels -> indeterminate. Unphased
       data consistent with several allele combinations (routine for CYP2D6 and
       NAT2) is genuine ambiguity: "positions present, allele unresolvable".
       Identical labels are not ambiguity, so they collapse.
    3. Either allele missing or unnamed -> indeterminate (structural no-call).
    4. Phenotype string says indeterminate -> indeterminate.
    5. Otherwise called.

    The diplotype *label* is never pattern-matched; see _INDETERMINATE_PHENOTYPES.
    """
    if not isinstance(entry, dict):
        raise PharmcatError(
            f"PharmCAT output for {gene} is not an object: {entry!r}"
        )

    diplotypes = entry.get("diplotypes")
    if diplotypes is None:
        diplotypes = []
    if not isinstance(diplotypes, list):
        raise PharmcatError(
            f"PharmCAT output for {gene} has non-list diplotypes: {diplotypes!r}"
        )
    for candidate in diplotypes:
        if not isinstance(candidate, dict):
            raise PharmcatError(
                f"PharmCAT output for {gene} has a non-object diplotype: "
                f"{candidate!r}"
            )

    indeterminate = GeneCall(gene, None, None, INDETERMINATE)
    if not diplotypes:
        return indeterminate

    labels = {candidate.get("label") for candidate in diplotypes}
    if len(labels) > 1:
        return indeterminate

    first = diplotypes[0]
    label = first.get("label")
    if label is not None and not isinstance(label, str):
        raise PharmcatError(
            f"PharmCAT output for {gene} has a non-string label: {label!r}"
        )

    if not (
        _allele_is_called(first.get("allele1"), gene)
        and _allele_is_called(first.get("allele2"), gene)
    ):
        return indeterminate

    phenotypes = first.get("phenotypes")
    if phenotypes is None:
        phenotypes = []
    if not isinstance(phenotypes, list):
        raise PharmcatError(
            f"PharmCAT output for {gene} has non-list phenotypes: {phenotypes!r}"
        )
    phenotype = phenotypes[0] if phenotypes else None
    if phenotype is not None and not isinstance(phenotype, str):
        raise PharmcatError(
            f"PharmCAT output for {gene} has a non-string phenotype: {phenotype!r}"
        )

    if _phenotype_is_indeterminate(phenotype):
        return indeterminate
    if not label:
        return indeterminate

    return GeneCall(gene, label, phenotype, CALLED)


def parse_phenotype_json(
    path: Path,
    uncovered_genes: Iterable[str],
    partially_covered_genes: Iterable[str],
) -> list[GeneCall]:
    """Translate PharmCAT phenotype output into GeneCall records.

    Coverage beats PharmCAT, always, and the precedence is strict:
    uncovered_genes > partially_covered_genes > PharmCAT output.

    uncovered_genes (Task 9: CoverageReport.genes_fully_uncovered) is a gene the
    array never informed at all -- no PharmCAT output about it can be treated as
    a call. partially_covered_genes (CoverageReport.genes_partially_covered) is
    worse than it looks: PharmCAT assumes reference at unobserved positions, so
    a gene with 1 of 40 positions covered still yields a confident "*1/*1
    Normal Metabolizer". That is the most dangerous output this software can
    produce, so partial coverage is reported as indeterminate.

    Both arguments accept any iterable of gene names and are normalized to
    frozenset. Neither has a default: a coverage guard that can be omitted is a
    guard that will be omitted.
    """
    uncovered = frozenset(uncovered_genes)
    partial = frozenset(partially_covered_genes) - uncovered

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError, OSError) as err:
        raise PharmcatError(f"could not parse PharmCAT output {path}: {err}") from err

    if not isinstance(payload, dict):
        raise PharmcatError(
            f"PharmCAT output {path} is not a JSON object "
            f"(got {type(payload).__name__})"
        )

    phenotypes = payload.get("phenotypes")
    if not isinstance(phenotypes, dict):
        raise PharmcatError(
            f"PharmCAT output {path} has no usable 'phenotypes' object "
            f"(got {type(phenotypes).__name__}); refusing to treat unparseable "
            f"output as an empty result"
        )

    calls: list[GeneCall] = []
    seen: set[str] = set()

    for gene, entry in phenotypes.items():
        seen.add(gene)
        if gene in uncovered:
            calls.append(GeneCall(gene, None, None, NOT_COVERED))
        elif gene in partial:
            calls.append(GeneCall(gene, None, None, INDETERMINATE))
        else:
            calls.append(_resolve(gene, entry))

    for gene in sorted(uncovered - seen):
        calls.append(GeneCall(gene, None, None, NOT_COVERED))
    for gene in sorted(partial - seen):
        calls.append(GeneCall(gene, None, None, INDETERMINATE))

    if not calls:
        raise PharmcatError(
            f"PharmCAT output {path} yielded no gene calls and no coverage "
            f"information was supplied; there is nothing to record"
        )

    return calls
