# Longitudinal Pharmacogenomic Record

**Date:** 2026-08-01
**Author:** Krrisha Patel
**Status:** Design approved

## Purpose

Pharmacogenomics has a delivery problem, not a science problem. Your genotype determines how you metabolize many common drugs, and CPIC publishes peer-reviewed, gene-drug-specific guidelines describing what to do about it. The evidence is settled and public.

What fails is persistence. A person is genotyped once. The result lands in a chart as a PDF. Three years later a different prescriber writes a drug that the genotype contraindicates, and nothing surfaces the conflict. The genotype is permanent; the drug exposures keep arriving for the rest of the person's life.

Existing tooling is batch-oriented. PharmCAT — maintained by the CPIC/PharmGKB group — takes a VCF, calls star alleles, and emits a report. It is excellent and this project does not attempt to replace it. But it is stateless: no memory of prior genotyping, no notion of time, no mechanism to revisit a stored result when a new drug is prescribed or when a guideline is revised.

This project builds that missing layer: a persistent, versioned pharmacogenomic record that is written once and re-evaluated on every subsequent drug query and every guideline update.

### What this is not

This is research and reference tooling. It makes no clinical claim about any specific person's care, and it is not clinical decision support. Two consequences shape the whole design:

- Output is framed as "the CPIC guideline for this gene-drug pair says X, here is the citation," never "this person should take dose Y."
- Guideline content is referenced by identifier and link, never redistributed as stored prose.

The second point is both a legal precaution and better engineering. Guidelines are revised; a stored copy is a stale copy.

## Scope of v1

Ingest a consumer raw-genotype file, produce a stored, versioned, citation-linked pharmacogenomic record, and answer drug queries against it.

**In scope**
- 23andMe / AncestryDNA raw genotype ingest, normalized to VCF
- Star-allele calling and phenotype assignment, delegated entirely to PharmCAT
- An immutable, append-only record store stamped with tool and guideline versions
- Drug query: given a stored record and a drug, return the applicable CPIC guidance reference
- Guideline-version diffing: when the pinned guideline version changes, report which stored records are affected

**Out of scope for v1**
- Any user-facing clinical recommendation or dose figure
- Clinical VCF ingest from sequencing labs (format variance is large; consumer arrays first)
- Multi-user auth, EHR/FHIR integration, prescriber workflows
- Structural variants and copy-number calling for `CYP2D6` (see Known Limitations)

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  Ingest                                                      │
│  23andMe / AncestryDNA raw text  ──▶  normalized VCF         │
│  (build detection, ref allele resolution, strand handling)   │
└──────────────────────────┬──────────────────────────────────┘
                           │ VCF
┌──────────────────────────┴──────────────────────────────────┐
│  Call  (thin wrapper — no calling logic of our own)          │
│  PharmCAT in Docker, pinned version, invoked as subprocess   │
│  ──▶ diplotypes + phenotypes + CPIC gene-drug pair IDs       │
└──────────────────────────┬──────────────────────────────────┘
                           │ parsed result
┌──────────────────────────┴──────────────────────────────────┐
│  Record  (the product)                                       │
│  Append-only store. Every entry stamped with:                │
│    subject id │ diplotype │ phenotype │ coverage status      │
│    pharmcat_version │ guideline_version │ ingested_at        │
└──────────────────────────┬──────────────────────────────────┘
                           │ stored phenotypes
┌──────────────────────────┴──────────────────────────────────┐
│  Re-evaluate                                                 │
│  drug query      ──▶ guidance reference + citation + link    │
│  guideline bump  ──▶ affected-record report                  │
└─────────────────────────────────────────────────────────────┘
```

The asymmetry is the point: genotype is ingested and called **once**; every later drug query is a cheap lookup against stored phenotypes, never a pipeline re-run.

## Component 1: Ingest

Consumer genotype files are messy in specific, known ways, and this is where most defects will live. Isolating it keeps that mess out of everything downstream.

Responsibilities:
- Detect vendor and genome build from the header; refuse rather than guess when ambiguous
- Map rsIDs to positions, resolve reference alleles, handle strand orientation
- Emit a valid VCF acceptable to PharmCAT

Interface: `path to raw file` → `path to VCF` + an ingest report listing every position that could not be confidently converted.

The ingest report is not diagnostic filler. Consumer arrays genotype a sparse subset of positions, and which positions are missing directly determines which alleles are callable. That must be visible downstream, not swallowed.

## Component 2: Call

A subprocess wrapper around PharmCAT, and nothing more. We never implement or "improve" star-allele calling.

- PharmCAT runs from a **pinned** container tag. Guideline and tool updates are deliberate, reviewed acts, never implicit.
- The wrapper parses PharmCAT's JSON output into our internal types and records the exact PharmCAT version used.
- If PharmCAT fails or emits anything unparseable, the record is **not** written. A partial record is worse than no record.

PharmCAT is MPL-2.0, which permits commercial and proprietary use around it; modifications to PharmCAT's own files must be published. We invoke it as a container and modify nothing, so this stays clean.

## Component 3: Record

Append-only. Records are never updated in place and never deleted.

Every entry carries, per gene:
- diplotype and assigned phenotype
- **coverage status** — one of `called`, `not_covered` (positions absent from the input), or `indeterminate` (positions present, allele unresolvable)
- `pharmcat_version`, `guideline_version`, `ingested_at`

Immutability is a correctness requirement, not a preference. A phenotype call is only meaningful relative to the tool and guideline versions that produced it. Overwriting a call destroys the ability to explain why the system once said something different — which is precisely the question that matters when guidance changes.

Re-ingesting the same subject appends a new entry. History is preserved and diffable.

## Component 4: Re-evaluate

Two directions:

**Drug query.** Given a subject and a drug, return the relevant gene(s), the stored phenotype, the CPIC gene-drug pair identifier, and a link to the current guideline. It returns a *reference*, not a recommendation.

**Guideline bump.** When the pinned guideline version advances, report which stored records involve gene-drug pairs whose guidance changed. This is the compounding behavior a stateless tool cannot offer: existing records gain value when guidelines are revised.

## The invariant that must not break

**Absence of a guideline and absence of coverage are different states, and must never collapse into a single "no interaction" answer.**

Three distinct outcomes, always distinguishable in the response:

| Outcome | Meaning |
|---|---|
| `guidance_found` | Gene called, CPIC guidance exists for this pair |
| `no_guidance_for_pair` | Gene called confidently, CPIC publishes nothing for this drug |
| `cannot_assess` | Gene `not_covered` or `indeterminate` — we do not know |

`cannot_assess` must never be rendered in a way that resembles reassurance. A consumer array that does not cover the relevant positions produces *no information*, and presenting that as "no interaction found" is the single most dangerous thing this system could do. This distinction is a test-suite obligation, not a UI nicety.

## Error handling

- **Ingest** — unknown vendor or ambiguous build: reject with a specific message. Never infer.
- **Call** — PharmCAT nonzero exit, timeout, or unparseable output: no record written, error surfaced with PharmCAT's stderr attached.
- **Record** — writes are atomic; a failed write leaves no partial entry.
- **Query** — unknown drug returns `no_guidance_for_pair` explicitly, never an empty success.

General rule: fail loudly and specifically. Silence is the failure mode with real-world consequences here.

## Testing

- **Ingest**: fixture files per vendor and build, including deliberately truncated and malformed inputs. Assert that unconvertible positions appear in the ingest report.
- **Call**: PharmCAT's own test VCFs, asserting our parser preserves its calls exactly. A golden-file test pinned to the PharmCAT version catches output-format drift on upgrade.
- **Record**: immutability properties — re-ingest appends rather than mutates; version stamps always present.
- **Re-evaluate**: one explicit test per invariant row above. The `cannot_assess` case gets a dedicated test asserting the response cannot be mistaken for "no interaction," including a case where the gene is entirely absent from a consumer array.
- **Guideline diffing**: two pinned guideline versions with a known delta; assert exactly the affected records are reported.

## Known limitations

Stated plainly because overstating capability is the main risk in this domain.

- **Consumer arrays are sparse.** They genotype a subset of positions. Many clinically relevant alleles are simply not callable from 23andMe data, and `CYP2D6` — among the most important PGx genes — depends on copy-number and structural variation that arrays cannot resolve. v1 reports these as `not_covered` rather than attempting inference.
- **No clinical validity.** Consumer genotype data is not clinically confirmed. This is a reference tool over self-supplied data.
- **Guideline lag.** Guidance is pinned; the tool is only as current as its pinned version. The guideline-bump report exists to make that lag visible instead of invisible.
- **Population bias.** Star-allele frequency and guideline evidence are unevenly distributed across ancestries. Guidance quality is not uniform across users.

## Open risks

1. **Guideline data licensing.** Verified this session: PharmCAT is MPL-2.0 and safe to build around. CPIC/PharmGKB data-use terms for commercial redistribution were **not** verifiable — `cpicpgx.org/license/` redirects to ClinPGx and the policy text did not render. v1 sidesteps this entirely by linking rather than embedding guideline prose. **Any future move to embed guideline content, or to commercialize, requires resolving this first** — one email to CPIC/PharmGKB, and counsel if the answer is not clearly permissive.
2. **No clinical design partner.** No pharmacist or clinical geneticist is reviewing this. It is the reason v1 is scoped as reference tooling with no clinical claims. Commercializing this without a clinician partner would be building clinically plausible software nobody clinically wants.
3. **PharmCAT output-format drift.** Version upgrades may change JSON structure. Mitigated by pinning plus golden-file tests.

## Future direction

Deliberately deferred, listed so the v1 boundary is legible:

- Clinical VCF ingest from sequencing labs
- FHIR / EHR integration
- Prescriber-facing surfaces — requires the clinician partner and a regulatory review of the CDS exemption
- `CYP2D6` structural-variant support via an additional caller

## Provenance

Design developed in conversation with Claude Opus 5 (Anthropic), which surveyed the PGx tooling landscape, verified PharmCAT's license and release cadence, identified that PharmCAT already covers the calling pipeline, and drafted this document.
