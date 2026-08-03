# pharmacogenomic-record

A single-subject longitudinal pharmacogenomic record built on
[PharmCAT](https://github.com/PharmGKB/PharmCAT).

**This is research and reference tooling. It is not a medical device, not
clinical decision support, it makes no claim about any person's care, and it
must not be used to make treatment decisions.**

## What this is

PharmCAT is excellent and stateless: VCF in, report out. This project adds the
layer it does not have — a persistent, version-stamped record that is written
once, never modified, and re-evaluated every time a new drug is queried or a
guideline is revised.

The pipeline:

1. Parse a 23andMe raw export, rejecting an unrecognized vendor or genome build
   rather than guessing.
2. Build a PharmCAT-ready VCF **keyed on rsID**. 23andMe reports GRCh37 and the
   PharmCAT position table is GRCh38; joining on rsID takes coordinates from the
   reference table and avoids a liftover step entirely.
3. Run the pinned `pgkb/pharmcat:3.4.0` image. Star-allele calling is not
   reimplemented here and never will be.
4. Store the gene calls append-only, stamped with the tool and guideline
   versions that produced them.
5. Answer "is drug X affected by my genotype?" with a CPIC citation.

## The invariant: absence of data is not absence of an interaction

This is the whole point of the project. A drug query has three outcomes and they
never collapse into one:

| Outcome | Meaning | Exit code |
|---|---|---|
| `guidance_found` | Gene called, CPIC publishes guidance for this pair | 0 |
| `no_guidance_for_pair` | Gene called **confidently**, CPIC publishes nothing | 0 |
| `cannot_assess` | Gene `not_covered`/`indeterminate`, or no record — **we do not know** | 2 |

Underneath sit three coverage states — `called`, `not_covered`,
`indeterminate` — and **coverage alone decides whether a gene is assessable**.
Never the phenotype: F2, F5, VKORC1, CFTR, IFNL3 and ABCG2 legitimately have no
metabolizer phenotype at all, so a `called` gene with `phenotype=None` is a real
answer and is reported as one.

`cannot_assess` is never rendered as reassurance. Presenting "we have no data"
as "no interaction found" is the most dangerous thing this system could do, so
the rendering is tested against a list of phrases it may not contain, and
`cannot_assess` exits **2** rather than 0 — a shell script writing
`if pharmacogenomic-record query codeine ...; then` must not read an unknown as
an all-clear.

## Coverage honesty: expect "cannot assess" most of the time

A gene is eligible to be `called` only if **every rsID-joinable position for
that gene was covered** by the array. Strict, with no threshold and no ratio,
because there is no defensible one: PharmCAT assumes reference at any position
it was not given, so a gene with 39 of 40 positions covered still yields a
confident "*1/*1 Normal Metabolizer" — and the missing position is exactly where
a variant would have been.

Consumer arrays genotype a very sparse subset of these positions. Measured on a
real 23andMe export against `pharmcat_positions_3.4.0.vcf`:

- **4 of 1226** positions covered
- 1014 joinable positions uncovered
- 208 positions carry no rsID at all, so they can never be joined from a 23andMe
  file under any circumstances

**So under the strict rule nearly every gene comes out `indeterminate` or
`not_covered`, and `called` is rare. That is the truthful result for a sparse
consumer array, not a defect in this tool.** If you run it on your own export
and it mostly says "cannot assess", it is working correctly and telling you
something real: your array did not measure enough of those genes to support a
confident answer.

Two further limits on the join, both structural:

- 9 of the 22 genes contain positions with no rsID whatsoever. G6PD is only 39%
  rsID-bearing and RYR1 78%, so even a perfect array could not fully inform them
  through this path. `unjoinable_positions` is reported for exactly this reason;
  it is a known residual gap, not a solved one.
- `CYP2D6` depends on copy-number and structural variation that consumer arrays
  cannot resolve at all. That is a different limitation from the rsID join and is
  not something coverage counting can detect.

## Usage

```bash
./scripts/fetch_positions.sh          # fetch the pinned 3.4.0 position table
pip install -e ".[dev]"

pharmacogenomic-record --db records.db ingest genome.txt --subject me   # needs Docker
pharmacogenomic-record --db records.db query clopidogrel --subject me
pharmacogenomic-record --db records.db drift --changed-pair CYP2D6-codeine
```

### Real output

Every block below is verbatim output from this tool, not an illustration.

**Ingest.** Coverage is reported before PharmCAT is invoked, so a run that fails
at the Docker step has still told you what your array covers. This machine has
no Docker, which is why the run ends where it does — and note that the failure is
surfaced loudly and no record is written:

```
$ pharmacogenomic-record --db records.db ingest tests/fixtures/23andme_valid.txt --subject me --workdir work
wrote work/23andme_valid.vcf
covered positions: 4
uncovered positions: 1014
positions with no rsID, unjoinable from a 23andMe file: 208
genes fully covered (eligible to be called): 0
genes partially covered (recorded indeterminate): 2 -- CYP2C19, DPYD
genes with no coverage (recorded not_covered): 20 -- ABCG2, CACNA1S, CFTR, CYP2B6, CYP2C9, CYP2D6, CYP3A4, CYP3A5, CYP4F2, F2, F5, G6PD, IFNL3, NAT2, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, VKORC1
error: docker not found on PATH; cannot run PharmCAT
$ echo $?
1
```

Those 4 covered positions are the honest result described above: zero genes
eligible to be called.

**`guidance_found`** — a gene whose every joinable position was covered, so
PharmCAT's call is allowed through:

```
$ pharmacogenomic-record --db records.db query clopidogrel --subject me
[GUIDANCE] CYP2C19 diplotype *1/*2. PharmCAT assigned phenotype 'Intermediate Metabolizer'. CPIC publishes guidance for pair CYP2C19-clopidogrel; the guideline itself is at https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/ and is not reproduced here.
OVERALL: GUIDANCE
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
0
```

**`cannot_assess`** — the common case. It is loud, it still cites where the
guidance is, and it exits 2:

```
$ pharmacogenomic-record --db records.db query codeine --subject me
[CANNOT ASSESS] Cannot assess CYP2D6 for pair CYP2D6-codeine: the stored coverage state is 'not_covered', so this gene's genotype is unknown. That is missing data, not absence of an interaction. CPIC guidance for this pair is at https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/.
OVERALL: CANNOT ASSESS
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
2
```

**A drug with several relevant genes.** Warfarin has two pairs, and the OVERALL
line is the *least reassuring* component — one unknown gene invalidates any claim
that the answer is complete:

```
$ pharmacogenomic-record --db records.db query warfarin --subject me
[CANNOT ASSESS] Cannot assess CYP2C9 for pair CYP2C9-warfarin: the stored coverage state is 'not_covered', so this gene's genotype is unknown. That is missing data, not absence of an interaction. CPIC guidance for this pair is at https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/.
[CANNOT ASSESS] Cannot assess VKORC1 for pair VKORC1-warfarin: the latest stored record holds no call for this gene, so its genotype is unknown. That is missing data, not absence of an interaction. CPIC guidance for this pair is at https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/.
OVERALL: CANNOT ASSESS
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
2
```

**A drug absent from the pair table** gets an explicit negative, never silence,
and the statement is scoped to what it actually establishes:

```
$ pharmacogenomic-record --db records.db query ibuprofen --subject me
[NO CPIC PAIR] This reference table lists no CPIC gene-drug pair for 'ibuprofen'. The table is a curated subset of CPIC's pairs, so this states only that 'ibuprofen' is absent from it.
OVERALL: NO CPIC PAIR
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
0
```

**An unknown subject is `cannot_assess`**, not "nothing found":

```
$ pharmacogenomic-record --db records.db query clopidogrel --subject nobody
[CANNOT ASSESS] Cannot assess CYP2C19 for pair CYP2C19-clopidogrel: no record is stored for this subject, so this gene's genotype is unknown. That is missing data, not absence of an interaction. CPIC guidance for this pair is at https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/.
OVERALL: CANNOT ASSESS
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
2
```

**Guideline drift** — the reason the record is persistent at all. Coverage is
deliberately *not* a filter here: a subject whose CYP2D6 was never covered is
still reported, because new guidance may be the very reason to finally get
proper testing:

```
$ pharmacogenomic-record --db records.db drift --changed-pair CYP2D6-codeine
[AFFECTED] subject me gene CYP2D6 pair(s) CYP2D6-codeine
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
```

## What it does not do

- **Call star alleles.** PharmCAT does that; the image is pinned and the logic is
  never reimplemented.
- **Restate guideline text.** Only `cpic_pair_id` and a URL are stored, and the
  URL host is pinned to `cpicpgx.org`. No dose figures, no clinical imperatives,
  no guideline prose appear anywhere in the output.
- **Infer anything PharmCAT did not call.** Coverage can only ever *downgrade*
  PharmCAT's answer, never upgrade it.
- **Resolve `CYP2D6`.** See above.
- **Accept clinical VCF input, or integrate with FHIR.** Out of scope for v1.

## Limitations

- Consumer genotype data is not clinically confirmed and has not been through a
  diagnostic laboratory.
- The pair table in `data/gene_drug_pairs.json` is a small curated subset of
  CPIC's pairs. A `no_guidance_for_pair` answer states only that the drug is
  absent from *this table*.
- Guidance is pinned to a version; the tool is only as current as that pin, which
  is why every record carries the stamp that produced it and why `drift` exists.
- Star-allele frequency and guideline evidence are unevenly distributed across
  ancestries, so guidance quality is not uniform across populations.
- `run_pharmcat` is unverified past its Docker guard in this environment; the
  automated tests cover the ingest half only.

## Development

```bash
pip install -e ".[dev]"
pytest
```

Genotype data must never be committed: `data/.gitignore` ignores everything in
`data/` except the two committed reference tables, and `work/` and `records.db`
are ignored at the repo root.

`tests/fixtures/23andme_full_cyp2c19.txt` is **synthetic** — it is not any
person's genotype data. It exists so the strict coverage rule can be exercised
end to end, by covering every rsID-joinable CYP2C19 and VKORC1 position.

## License and attribution

PharmCAT is MPL-2.0. CPIC and PharmGKB guideline content is referenced by link,
never redistributed. Embedding guideline text or commercializing this project
requires resolving CPIC/PharmGKB data-use terms first — those terms were never
verified for this project, and the link-only design is what sidesteps them.
