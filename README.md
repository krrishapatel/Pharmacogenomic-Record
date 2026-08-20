# Pharmacogenomic-Record

[![CI](https://github.com/krrishapatel/Pharmacogenomic-Record/actions/workflows/ci.yml/badge.svg)](https://github.com/krrishapatel/Pharmacogenomic-Record/actions/workflows/ci.yml)

A single-subject longitudinal pharmacogenomic record built on
[PharmCAT](https://github.com/PharmGKB/PharmCAT).

**This is research and reference tooling. It is not a medical device, not
clinical decision support, it makes no claim about any person's care, and it
must not be used to make treatment decisions.**

## What this is

PharmCAT is excellent and stateless: VCF in, report out. This project adds the
layer it does not have: a persistent, version-stamped record that is written
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
| `cannot_assess` | Gene `not_covered`/`indeterminate`, or no record, so **we do not know** | 2 |

Underneath sit three coverage states: `called`, `not_covered`, and
`indeterminate`. **Coverage alone decides whether a gene is assessable**.
Never the phenotype: F2, F5, VKORC1, CFTR, IFNL3 and ABCG2 legitimately have no
metabolizer phenotype at all, so a `called` gene with `phenotype=None` is a real
answer and is reported as one.

`cannot_assess` is never rendered as reassurance. Presenting "we have no data"
as "no interaction found" is the most dangerous thing this system could do, so
the rendering is tested against a list of phrases it may not contain, and
`cannot_assess` exits **2** rather than 0. A shell script writing
`if pharmacogenomic-record query codeine ...; then` must not read an unknown as
an all-clear.

## Coverage honesty: expect "cannot assess" most of the time

A gene is eligible to be `called` only if **every rsID-joinable position for
that gene was covered** by the array. Strict, with no threshold and no ratio,
because there is no defensible one: PharmCAT assumes reference at any position
it was not given, so a gene with 39 of 40 positions covered still yields a
confident "*1/*1 Normal Metabolizer", and the missing position is exactly where
a variant would have been.

Two facts about `pharmcat_positions_3.4.0.vcf` itself are fixed and verified
(see `tests/test_reference_data.py`), independent of any array:

- the table has **1226 positions**;
- **208 of them carry no rsID at all**, so they can never be joined from a
  23andMe file under any circumstances, no matter how good the array is.

A consumer array genotypes only a very sparse subset of the 1018 rsID-bearing
positions. So under the strict rule nearly every gene comes out `indeterminate`
or `not_covered`, and `called` is rare. That is the truthful result for a
sparse consumer array, not a defect in this tool. If you run it on your own
export and it mostly says "cannot assess", it is working correctly and telling
you something real: your array did not measure enough of those genes to support
a confident answer. This project does not ship anyone's real genome, so the
concrete covered/uncovered counts shown below are from the small bundled test
fixtures, not a measurement of any personal array; a real export would cover
more positions than a fixture, but still a sparse subset of the 1018.

Two further limits on the join, both structural:

- 9 of the 22 genes contain positions with no rsID whatsoever. G6PD is only 39%
  rsID-bearing and RYR1 78%, so even a perfect array could not fully inform them
  through this path. `unjoinable_positions` is reported for exactly this reason;
  it is a known residual gap, not a solved one.
- Hemizygous calls are not encoded, which for a male sample means **no G6PD
  position is usable at all**. G6PD is the only X-linked gene in the table, and
  23andMe reports one allele per position for male non-PAR chrX, so all 67 of its
  rsID-joinable positions arrive as single letters and none becomes a VCF row.
  This is a limitation of this pipeline and not of VCF: a haploid `GT` of `0` or
  `1` is valid VCF 4.2 and is what those calls should become. It is left
  unimplemented because whether the pinned PharmCAT image accepts a haploid GT
  for G6PD has not been verified here, and feeding it something unverified is
  worse than reporting the gap. What was fixed is the reporting: those positions
  used to be dropped during parsing, so they were indistinguishable from
  positions the array never measured and G6PD came back `not_covered`, meaning
  "your array said nothing about this gene" when it had in fact measured every
  position. They are now carried through to `hemizygous_rsids` and printed as
  `positions measured as a single allele, not encodable here`. They still count
  as uncovered, so no call is affected; the difference is that the reason is now
  the true one.
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

### Example output

Every block below is verbatim output from this tool, not an illustration. The
inputs are the small **bundled test fixtures** under `tests/fixtures/`, not
anyone's real genome. The command line names the fixture in each case; they are
not a measurement of any personal array.

The two **ingest** blocks are reproducible on any checkout: their coverage
numbers are properties of the named fixtures and are re-derived by the test
suite, so they do not need Docker. The **`query`** and **`drift`** blocks that
follow read from a *stored record*, and writing one requires `run_pharmcat`
(Docker), which is not available here. They were produced on a Docker-capable
machine after ingesting the synthetic `23andme_full_cyp2c19.txt` fixture end to
end; on a fresh checkout without that stored record, `query clopidogrel` instead
prints a `CANNOT ASSESS` line ("no record is stored for this subject"). They are
shown here to illustrate the three query outcomes, not as output you can
regenerate without Docker.

**Ingest.** Coverage is reported before PharmCAT is invoked, so a run that fails
at the Docker step has still told you what your array covers. This machine has
no Docker, which is why the run ends where it does. Note that the failure is
surfaced loudly and no record is written. The `23andme_valid.txt` fixture is a
10-line sample that covers just 4 positions:

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

Those 4 covered positions (out of the 1226 in the table, 208 of which carry no
rsID) leave zero genes eligible to be called, the sparse-array result the
strict rule is built for. A larger fixture,
`tests/fixtures/23andme_full_cyp2c19.txt`, which is **synthetic**, engineered
to cover every rsID-joinable CYP2C19 and VKORC1 position, reports 36 covered
positions and two fully covered genes, showing the other side of the rule:

```
$ pharmacogenomic-record --db records.db ingest tests/fixtures/23andme_full_cyp2c19.txt --subject me --workdir work
wrote work/23andme_full_cyp2c19.vcf
covered positions: 36
uncovered positions: 982
positions with no rsID, unjoinable from a 23andMe file: 208
genes fully covered (eligible to be called): 2 -- CYP2C19, VKORC1
genes partially covered (recorded indeterminate): 0
genes with no coverage (recorded not_covered): 20 -- ABCG2, CACNA1S, CFTR, CYP2B6, CYP2C9, CYP2D6, CYP3A4, CYP3A5, CYP4F2, DPYD, F2, F5, G6PD, IFNL3, NAT2, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1
error: docker not found on PATH; cannot run PharmCAT
$ echo $?
1
```

**`guidance_found`**: a gene whose every joinable position was covered, so
PharmCAT's call is allowed through:

```
$ pharmacogenomic-record --db records.db query clopidogrel --subject me
[GUIDANCE] CYP2C19 diplotype *1/*2. PharmCAT assigned phenotype 'Intermediate Metabolizer'. CPIC publishes guidance for pair CYP2C19-clopidogrel; the guideline itself is at https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/ and is not reproduced here.
OVERALL: GUIDANCE
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
$ echo $?
0
```

**`cannot_assess`**: the common case. It is loud, it still cites where the
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
line is the *least reassuring* component. One unknown gene invalidates any claim
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

**Guideline drift**: the reason the record is persistent at all. Coverage is
deliberately *not* a filter here: a subject whose CYP2D6 was never covered is
still reported, because new guidance may be the very reason to finally get
proper testing:

```
$ pharmacogenomic-record --db records.db drift --changed-pair CYP2D6-codeine
[AFFECTED] subject me gene CYP2D6 pair(s) CYP2D6-codeine
This output is a reference to published CPIC citations, keyed on stored gene calls. It is not a medical device, not clinical decision support, and not a basis for any treatment decision.
```

A mistyped pair id gets a **warning**, never the reassuring "touches nobody"
negative a real-but-unmatched pair would get. An id that matched nothing was
never checked, and that is not the same as a revision touching nobody:

```
$ pharmacogenomic-record --db records.db drift --changed-pair CYP2D6-codiene
WARNING: changed pair(s) CYP2D6-codiene match no CPIC gene-drug pair in this table. This is almost certainly a mistyped id; nothing was checked for it, which is not the same as its revision touching nobody. Verify the id against cpicpgx.org.
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
- The pair table in `src/pharmacogenomic_record/data/gene_drug_pairs.json` is a small curated subset of
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

Genotype data must never be committed: the reference tables ship as package
data under `src/pharmacogenomic_record/data/` (so they are found through
`importlib.resources` and travel inside the wheel, not via a path relative to
the repo root), and `src/pharmacogenomic_record/data/.gitignore` ignores
everything in that directory except the two committed reference tables. `work/`
and `records.db` are ignored at the repo root.

`tests/fixtures/23andme_full_cyp2c19.txt` is **synthetic**. It is not any
person's genotype data. It exists so the strict coverage rule can be exercised
end to end, by covering every rsID-joinable CYP2C19 and VKORC1 position.

## License and attribution

PharmCAT is MPL-2.0. CPIC and PharmGKB guideline content is referenced by link,
never redistributed. Embedding guideline text or commercializing this project
requires resolving CPIC/PharmGKB data-use terms first. Those terms were never
verified for this project, and the link-only design is what sidesteps them.
