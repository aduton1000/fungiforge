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
7 min 10 s, 0.4 CPU-h; triage on 200,000 of 11.2 M read pairs: 93.5 % bacterial, 3.8 % unclassified,
top taxon *Staphylococcus aureus* (17.3 %); master row present with `stages_failed = gate:skipped`,
`gate = skipped(read_triage:non_fungal)`; no assembly, annotation or downstream task ran. On the
production code the same isolate went through assembly and hours of annotation.

### 1.2 Multi-locus identification (W2.3)

Pending on the dev deployment (needs `fetch_references.sh markers mlst busco` on the data root):
CEA10 expected *A. fumigatus* at high confidence with at least one agreeing secondary locus, an
afumigatus MLST ST and eurotiales_odb10 completeness; DF-005 expected ITS *A. flavus* with CaM and
BenA reported as ties (*A. flavus*/*A. oryzae*) and TEF1/RPB2 deciding.

Real-tool check, 2026-09-17, base image on a workstation, CEA10 nuclear assembly from the v0.1.0
run, reference sets fetched the same day (CaM 3,829, BenA 8,031, LSU 55,721 type-material records;
PubMLST afumigatus and calbicans):

| Locus | Extracted | Best type-material hit |
|---|---|---|
| CaM | contig_45, 2 HSPs (intron), 1,021 bp | *A. fumigatus* KACC 41143 / NRRL 163, 100 % over 505 bp |
| BenA | contig_20, 3 HSPs, 1,843 bp | *A. fumigatus* CBS 133.61 / NRRL 163, 100 % over 390–444 bp |
| TEF1 / RPB2 | found (TEF1 1 HSP; RPB2 minus strand, 3.6 kb) | sets not fetched in this check |
| LSU D1/D2 | contig_8, 900 of 3,363 bp | 99.9 % to eight section *Fumigati* type strains: genus-level only |
| MLST | — | afumigatus ST5, 7/7 exact alleles |

Call without ITS: *A. fumigatus*, medium, `CaM+BenA(secondary-only)`, lineage eurotiales_odb10.

### 1.3 Annotation (W2.5)

Pending on the dev deployment: CEA10 and DF-005 rerun with `--genemark_dir` / `--genemark_key`
(GeneMark-ES joins EVM), species-aware training (eurotiomycetes set, *A. fumigatus* / *A. oryzae*
Augustus seeds) and eggNOG-mapper; record gene counts and the new coverage columns (`pct_pfam`,
`pct_go`, `pct_eggnog`) against the v0.1.0 numbers (9,630 and 12,473 proteins).

### 1.4 Extras (W2.6)

Pending on the dev deployment with the rebuilt image, the dbCAN/PHI-base/EffectorP data and the
site-built SignalP 6 image: CEA10 expected MAT1-1 (published for CEA10) and `haploid_like`;
DF-005 haploid; secretome, CAZyme and effector counts within the ranges published for the two
species.

### 1.5 Mitochondrial genome (W2.7)

Pending on the dev deployment: CEA10 expected one circular contig of about 31 kb (Af293:
30,696 bp) with 15/15 core genes and a mito/nuclear depth ratio in the tens; DF-005 (Illumina)
mitogenome separated and annotated. In-image check on the local CEA10 assembly: 15/15 core genes,
cox1 intron 2,011 bp (Af293 2,020 bp), rRNAs 100 %, redundant Flye fragments recognised (size
estimate 29.4 kb), copy ratio 32 (details in UPGRADE_PLAN W2.7).

### 1.6 Mobile elements (W2.8)

Pending on the dev deployment after `fetch_references.sh genomad`: `te_percent` populated for CEA10
(about 4.7 % from the earlier RepeatMasker table) and DF-005, geNomad status `ok`, RVDB screen
categories on both isolates. In-image RVDB screen on the local CEA10 assembly: 30 loci, one
partitivirus coat-protein EVE candidate (81 % identity), te_percent 4.67 from the RepeatMasker table
(details in UPGRADE_PLAN W2.8).

### 1.7 Novelty (W2.9)

Pending on the dev deployment after `fetch_references.sh genomes`: CEA10 and DF-005 `known_species`
with their nearest reference genomes (Af293 / NRRL 3357) above 98 % ANI. In-image skani check on
the local CEA10 assembly: 99.59 % ANI to Af293 (aligned fraction 94 %), known_species.

### 1.8 Cohort phylogenomics and clonality (W3.1)

Pending on the dev deployment: the study batch's species clusters against the identifications,
the supermatrix tree, and SNP distances within the *A. flavus* cluster; the published CEA10 / C87 /
C6 / E142 relationships as the topology check. In-image mini cohort (CEA10, Af293, NRRL 3357): CEA10 and
Af293 cluster at 99.59 % ANI with 42,457 core SNVs over 27.4 Mb; NRRL 3357 separate (UPGRADE_PLAN W3.1).

### 1.9 BGC families and mycotoxin flags (W3.2)

Pending on the dev deployment (antiSMASH now runs KnownClusterBlast): DF-005 expected to flag the
aflatoxin cluster, CEA10 gliotoxin plus the other known *A. fumigatus* clusters; the batch's shared
gene-cluster families. In-image family clustering of the 40 CEA10 regions (duplicated as a second
isolate): 40 families, each shared by both copies, none merged (UPGRADE_PLAN W3.2).

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
pipeline outcome for the positives (W2.4 code): `cyp51A_TR = TR34` with the read check
`agrees_with_assembly`, `L98H` known mutation with `read_support.agrees = true`, azole in
`resistant_classes`, confidence `high` (read-confirmed even without Illumina); for 157DB3: no TR
(reads `agrees_with_assembly` at the wild-type site), no L98H, the five background substitutions
reported only where they sit on a panel position (F46Y, M172V, N248T, D255E, E427K are all FungAMR
positions and will appear as `associated_unvalidated` or wild-type-agreeing calls), so
`resistant_classes = none`; `resistance_read_support` lists the confirmed positions.

In-image check of the read-level code (2026-09-17, simulated 250-bp reads on the real cyp51A
locus): C87 published assembly → L98H confirmed by 25 reads (96 %), TR34 agrees; wild-type
assembly with reads carrying L98H and an extra TR34 unit → L98H reported as reads-only (97 %) and
TR34 as `reads_have_extra_copy` (12/12 spanning reads), i.e. the C87 lesson is now detectable.

*Candida* controls (*C. auris* FKS1 S639F, *C. albicans* ERG11) are deferred to the resistance
upgrade (W2.4), where the caller gains read-level genotyping for those genes; candidate sources:
Spruijtenburg et al. 2023 (*Mycoses*, PMID 37712885) for *C. auris*.

## 4. Change log

- 2026-09-17 — DF-003 (bacterial) stopped at read triage on the dev deployment in 7 minutes: the gate validated on real data.
- 2026-09-17 — first dev-deployment runs of CEA10 and DF-005 compared: all scientific checks within tolerance; validator stage-lookup bug found and fixed; L30 (purge_dups never ran) and L31 (NanoPlot failures) opened.
- 2026-09-16 — suite created: expected files for CEA10 and DF-005, validate/benchmark/control tools with unit and minimap2 integration tests, controls chosen and fetch script written; C87 lesson recorded.
