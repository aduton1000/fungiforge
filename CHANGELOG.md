# Changelog

All notable changes to FungiForge. Versions follow the `docs/UPGRADE_PLAN.md` work items; the
master-table columns are append-only, so a newer release never breaks a Layer-2 consumer.

## [0.2.0] — unreleased (branch `develop`)

The pipeline version, the Python package and the default image tag are all 0.2.0 from this point,
so a dev deployment builds `fungiforge-0.2.0.sif` and `bin/lock_env.sh` finds the matching tag.

The correctness and completeness release: every stage that silently did nothing now does what the
manual says, or records why it could not. Validated on *A. fumigatus* CEA10 (hybrid), *A. flavus*
DF-005 (Illumina-only) and a bacterial study isolate; see `docs/validation.md`.

### Added
- **Status contract and provenance** (W0.1, W0.2): every stage writes `status`, per-tool exit
  codes and versions; `DB_CHECK` fails a run in seconds when a database is missing;
  `provenance.json` records commit, parameters, container digests and every tool version.
- **Test framework and CI** (W0.3): unit tests over every helper, nf-test for the DAG and each
  module stub, `nextflow lint` and shellcheck in GitHub Actions.
- **Validation suite** (W1.1): `fungiforge validate` against expected-results files,
  `benchmark_assembly.py`, read-level control verification (`control_genotype.py`).
- **Gating** (W2.1): non-fungal or failed-QC isolates stop before the heavy stages with the reason
  on their master row; `--force_all` overrides.
- **Read QC** (W2.2): trimmed reads are used downstream, Kraken2 read triage before assembly, and a
  KMC + GenomeScope2 genome-size / heterozygosity / ploidy / coverage profile.
- **Identification** (W2.3): CaM, BenA, TEF1, RPB2 and LSU extracted from the assembly and searched
  against NCBI type-material sets; a species call needs ITS plus an agreeing secondary locus; ties
  and discordance are flagged; fungal PubMLST MLST; species-aware BUSCO lineage (stage 08b).
- **Resistance** (W2.4): panel merged with the FungAMR catalogue (evidence tiers), read-level
  genotyping of every hotspot (allele frequency, zygosity, alleles the assembly missed), the
  cyp51A promoter repeat re-checked in the reads, and locus/genome copy-number ratios.
- **Annotation** (W2.5): GeneMark-ES from a licensed directory, species-aware Augustus and BUSCO
  training, eggNOG-mapper (stage 07b) and InterProScan (stage 07c) as their own stages, and
  annotation-coverage statistics.
- **Auxiliary biology** (W2.6): SignalP 6 (site-built image) secretome, EffectorP 3, dbCAN CAZymes,
  PHI-base virulence, Pfam-based mating-type typing with GenBank synteny, nQuire ploidy.
- **Organelle** (W2.7): the mitochondrial genome is separated from the assembly and annotated
  (core genes, introns, rRNAs, circularity, copy number, heteroplasmy) in the new stage 04b.
- **Mobile elements** (W2.8): `te_percent` and the TE landscape from RepeatMasker, geNomad on the
  whole assembly, an RVDB mycovirus / endogenous-viral-element screen, mitochondrial homing
  endonucleases.
- **Novelty** (W2.9): a staged reference genome set, skani ANI against it, and verdicts combining
  ANI with ITS and the secondary loci.
- **Cohort analysis** (W3.1, W3.2): run-level stage 16 — ANI species clusters, a single-copy BUSCO
  supermatrix tree (IQ-TREE 3), assembly-based core SNP distances and clonal groups, and
  gene-cluster families (BiG-SCAPE 2, or shared-protein clustering as fallback).
- **Reporting** (W3.3): a JSON schema for the master table with a validator, a Jinja2 per-isolate
  report, and run-level stage 17 writing the merged `master_fungi.tsv`, `cohort_report.html`,
  `cohort_summary.json` and MultiQC custom content.
- **Pre-flight** (W4.1): `fungiforge check --samplesheet` validates ids, read pairing, file
  existence and metadata consistency before a run; duplicate ids now fail at launch.

### Fixed
- Identification: a locus call is anchored on the longest alignment among hits that reach the
  species threshold; a shorter type-material record of higher identity can widen it into a tie
  but no longer replaces it (DF-005: a 531-bp *A. kambarensis* CaM record had displaced the
  *A. flavus* neotype and the isolate was reported *Aspergillus* sp.). Retired names are folded
  into accepted species through `fungiforge/resources/species_synonyms.tsv` (cited rows).
- Resistance: a curated panel row that names a species applies to that species only (genus rows
  still cover the genus); every call records the reference used, and residue calls against a
  reference from another species are reported as a cross-species screen and never counted as
  resistance (DF-005, *A. flavus*, had four *A. fumigatus*-numbered interspecies differences
  reported as azole resistance). FungAMR evidence from another species of the genus is an
  association, not a known mutation.
- `hpc_install.sh` refuses to update a checkout with local modifications and no longer hides a
  failed fast-forward behind `|| true`; a branch ref now deploys exactly what origin has.
- Illumina-only (and ONT-only) isolates no longer report `polish:skipped` in `stages_failed`: the
  polish stage's passthrough is its job for those isolates, so its Polypolish skip is recorded as
  not applicable (`ff_skip --not-applicable`, listed under `not_applicable` in the stage JSON) and
  the status stays `ok`. A stage that ran no tool because a database was absent still reads
  `skipped`. A static test now rejects any module branch that plain-skips without running a tool.
- Haplotig purging ran in an image without `purge_dups`, so it exited 127 on every long-read
  isolate and left the assembly stage `partial`; it is now its own stage in the base image (L30).
- Stage 05's QUAST call could never succeed and made every run's QC `partial` (L29).
- NanoPlot exited 1 on real ONT reads with plotly ≥ 6 (L31).
- Medaka was given a model name that need not match the reads; it is now read from the read
  headers and checked against the installed models (L27).
- Every containerised task gets its own `/tmp`, after a node tmpfs overflow broke annotation and
  silently failed barrnap (L28).
- The mitochondrial split wrote an empty placeholder (L14); geNomad was never staged or run (L15);
  the reference genome set for novelty was never staged (L16).

## [0.1.3] — 2026-09-16
- Per-task `/tmp` bind for every containerised process (node tmpfs overflow).

## [0.1.2] — 2026-09-16
- Annotation task writes its scratch under the work directory.

## [0.1.1] — 2026-09-15
- Removed nf-test artefacts from the release tag; `.gitignore` hardened.

## [0.1.0] — 2026-09-11
First working pipeline: ONT / Illumina / hybrid routing, assembly, polishing, decontamination,
QC, repeat masking, Funannotate annotation, ITS identification, antifungal-resistance calling
with the cyp51A TR detector, BGC detection, novelty, per-isolate reports and the master table.
Validated end to end on *A. fumigatus* CEA10 and *A. flavus* DF-005.
