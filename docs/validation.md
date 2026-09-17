# Validation record

What has been demonstrated, on which data, with which version. Every entry here is backed by a
run whose provenance file exists; "designed" features are not listed as validated. Companion
tooling (W1.1 of `UPGRADE_PLAN.md`):

| Tool | Purpose |
|---|---|
| `bin/validate_run.py` (`fungiforge validate`) | compares a finished run with `test/expected/<sample>.json` (tolerances per metric) and prints a pass/fail table; exit 1 on any failure |
| `bin/benchmark_assembly.py` | assembly vs reference genome with minimap2: aligned fractions, NGA50, breakpoints, mismatches and indels per 100 kb, unaligned contigs |
| `bin/control_genotype.py` + `test/controls/` | verifies a control isolate's genotype **from its reads** at a wild-type locus (promoter insertion and point-mutation checks) before the control is trusted |
| `bin/fetch_controls.sh` | downloads the public controls, runs the read-level check, writes a samplesheet row per verified control |
| `bin/fetch_references.sh benchmarks` | reference genomes for the benchmarks (A1163 for CEA10, NRRL 3357 for *A. flavus*) |

## 1. Reference isolates (expected-results files)

| Isolate | Data | Path | Expected file | Recorded from |
|---|---|---|---|---|
| *A. fumigatus* CEA10 | ONT SRR28036069 + Illumina SRR13127874 | hybrid | `test/expected/AfumCEA10.json` | shared-cluster run, v0.1.0, 2026-09-11: 29.9 Mb, 43 contigs, BUSCO 99.7 %, 9,630 proteins, 50 BGC, ITS 100 % *A. fumigatus*, cyp51A wild-type, 15/15 stages |
| *A. flavus* DF-005 (study isolate) | Illumina 151 bp PE | short-read | `test/expected/ASSARM-PHI-DF-005.json` | shared-cluster run, v0.1.0, 2026-09-14: 37.8 Mb, 70 contigs, N50 1.04 Mb, BUSCO 99.7 %, 12,473 proteins, 102 BGC, 6 human contigs removed |

Tolerances are recorded inside each expected file. Run after any pipeline change:

```bash
fungiforge validate --results ~/runs/dev/cea10/results --sample AfumCEA10
fungiforge validate --results ~/runs/dev/df005/results --sample ASSARM-PHI-DF-005
```

Status (2026-09-17, first runs of the `develop` code on the dev deployment, checkout `cc534dc`,
i.e. W0.1–W0.3 + the /tmp fix, before W1.1/W2.x):

| Check | CEA10 (hybrid) | DF-005 (Illumina-only) |
|---|---|---|
| species / confidence | *A. fumigatus* / high, ITS 100 % | *A. flavus* (ITS) |
| assembly | 29,986,658 bp, 48 contigs, N50 2,646,760, BUSCO 99.7 %, qc_pass | 37,804,321 bp, 70 contigs, N50 1,038,584, BUSCO 99.7 %, qc_pass (identical to the recorded run) |
| proteins | 9,607 (expected 9,630 ± 5 %) | 12,464 (expected 12,473 ± 5 %) |
| BGC regions | 52 (expected 50 ± 15 %) | 102 (expected 102 ± 15 %) |
| resistance | none; 0 known mutations; cyp51A located, promoter wild-type | none |
| decontamination | verdict fungal | verdict fungal, 0.04 % of bases removed |
| novelty | known_species | — |
| provenance | 27 tool versions, images hashed; all stage JSONs carry status/versions | same |
| stages not `ok` | readqc:partial (NanoPlot exit 1, L31), assemble:partial (purge_dups not in the Flye image, L30), assembly_qc:partial (phantom QUAST, fixed) | assembly_qc:partial (phantom QUAST, fixed) |

Validator result with the fixed lookup: CEA10 20/21, DF-005 15/16 — the single failure on each is
`stages_all_ok`, i.e. the `partial` stages listed above, all from code that has since changed except
purge_dups (L30). Every scientific check passes within tolerance on both isolates. The first pass of the validator
reported the assembly-QC metrics as missing because it looked the stage JSON up by stage name
(`assembly_qc`) while the file is `<sample>.assemblyqc.json`; fixed (`6f8fd55`), a real bug the
first real run exposed. CEA10 wall time 2 h 09 m on the production run; on the dev install the
rerun after the node fixes reused 12 cached tasks. DF-005: 9 h 09 m, 215 CPU-h (annotation).

Known deviation for runs made before 2026-09-16: the assembly-QC stage reports `partial` because
its optional QUAST call could never run (QUAST is in none of the images); the call has been
removed, so `stages_failed` is `none` only for runs of the current code.

### 1.1 Gate on a bacterial study isolate (W2.1 / W2.2)

DF-003 (Illumina-only; the plate's read triage had flagged it bacterial) on the dev deployment with
the W2.2 code (`8b0571e`), 2026-09-17: read QC → read triage → gate → report → provenance, 6 tasks,
7 min 10 s, 0.4 CPU-h; master row present with `stages_failed = gate:skipped`,
`gate = skipped(read_triage:non_fungal)`; no assembly, annotation or downstream task ran. On the
production code the same isolate went through assembly and hours of annotation.

## 2. Reference-genome benchmarks

| Isolate | Reference | Status |
|---|---|---|
| CEA10 | *A. fumigatus* A1163 (GCA_000150145.1) | pending: `bin/benchmark_assembly.py` on the dev-deployment CEA10 assembly |
| DF-005 | *A. flavus* NRRL 3357 (GCA_009017415.1) | pending |

## 3. Resistance controls

### 3.1 The C87 lesson (2026-09-15)

Isolate C87 (Hemmings, Rhodes & Fisher 2023) is published as TR34/L98H, and its published assembly
(GCA_949125185) does carry TR34 and L98H: `bin/cyp51a_TR.py` reports TR34 with two copies on it
(regression fixture `test/fixtures/cyp51A/`). The reads deposited under the same name (ONT
ERR10820709, Illumina ERR9791656) do **not**: 0 of 84 ONT and 1 of 76 Illumina reads spanning the
site carry an insertion, both read sets agree with each other, and the genome assembled from them
differs from the published assembly by 0.64 % with a 34-bp gap at the TR site. They are a
different isolate (cyp51A Y121F, T289A, G448S; wild-type promoter). The pipeline called that
genotype correctly at high confidence. Consequence: **a control is only accepted after its
genotype has been confirmed from its own reads**, which is what `control_genotype.py` and
`fetch_controls.sh` enforce.

### 3.2 Controls chosen (A. fumigatus cyp51A)

From the 2026 poultry-farm study (PRJNA1494181; ONT PromethION, Dorado high-accuracy basecalling;
genotypes from the paper's Table 3):

| Control | Run | Claimed genotype | Read-level check | Pipeline run |
|---|---|---|---|---|
| afum_TR34_L98H_226BB2 | SRR39591703 | TR34, L98H (+ F46Y M172V N284T D255E E427K) | pending | pending |
| afum_TR34_L98H_160CME3 | SRR39591704 | TR34, L98H (same background) | pending | pending |
| afum_wt_157DB3 | SRR39591705 | no TR, no L98H, same background | pending | pending |
| afum_CEA10 | SRR28036069 + SRR13127874 | wild-type | pending | passed as reference isolate (assembly-level) |

The three poultry-farm isolates share the same cyp51A background substitutions, so the TR34 and
L98H calls are tested against a matched wild-type rather than an unrelated strain. Expected
pipeline outcome for the positives: `cyp51A_TR = TR34`, `L98H` known mutation, azole in
`resistant_classes`, confidence `provisional_ont_only` (no Illumina); for 157DB3: no TR, no L98H,
the five background substitutions reported as novel hotspot variants only where they sit on a
panel hotspot (none do), so `resistant_classes = none`.

*Candida* controls (*C. auris* FKS1 S639F, *C. albicans* ERG11) are deferred to the resistance
upgrade (W2.4), where the caller gains read-level genotyping for those genes; candidate sources:
Spruijtenburg et al. 2023 (*Mycoses*, PMID 37712885) for *C. auris*.

## 4. Change log

- 2026-09-17 — DF-003 (bacterial) stopped at read triage on the dev deployment in 7 minutes: the gate validated on real data.
- 2026-09-17 — first dev-deployment runs of CEA10 and DF-005 compared: all scientific checks within tolerance; validator stage-lookup bug found and fixed; L30 (purge_dups never ran) and L31 (NanoPlot failures) opened.
- 2026-09-16 — suite created: expected files for CEA10 and DF-005, validate/benchmark/control tools with unit and minimap2 integration tests, controls chosen and fetch script written; C87 lesson recorded.
