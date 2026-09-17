# FungiForge v0.2 — upgrade plan and progress record

**Living document.** This is the single source of truth for the v0.2 upgrade: what is being built, why, how each piece is validated, and where the work stands. It is updated every time an item changes state, and it is the first thing read when work resumes. Status board first, detail below.

Started 15 September 2026 from v0.1.0 (commit `92a484e`, tagged `v0.1.0`). Development happens on branch `develop`; `main` stays at v0.1.0 until the cluster's production runs (77-isolate ASSARM-PHI batch, C87 positive control) are complete, then v0.2.0 is merged, tagged and deployed with `bin/hpc_install.sh --ref v0.2.0 --rebuild-images`.

## Ground rules

1. **Every item ships four things**: the implementation, an automated test (unit test for Python, stub or nf-test for workflow logic), documentation (manual + README where user-facing), and a validation run on real data with recorded results. An item is not done until all four exist.
2. **No silent failure, anywhere.** Every stage JSON carries `status` (`ok | partial | failed | skipped`) and the exit status of each tool it ran. A tool failure is either fatal to the task or recorded as `failed` with the log tail; `|| true` is never used to make a stage look successful.
3. **Reproducible by construction.** Versioned image tags (no `latest`), exact conda pins, a per-run `provenance.json` listing every tool version, image digest and database release.
4. **Honest documentation.** The manual describes what the code does today. Aspirational text is labelled as roadmap.
5. **Validation data are fixed**: CEA10 hybrid (ONT SRR28036069 + Illumina SRR13127874), ASSARM-PHI-DF-005 Illumina-only (*A. flavus*), C87 hybrid TR34/L98H positive (ERR10820709 + ERR9791656), plus the controls added in W1.1. Expected results live in `test/expected/` and `bin/validate_run.py` compares a run against them with tolerances.
6. **Backwards compatibility of the master table**: columns are append-only; the Layer-2 scripts are updated in the same change as any new column.
7. Work in small commits on `develop`, each referencing the item ID (`W2.4:`), with the status board updated in the same commit.

## Status board

Legend: `todo` · `in progress` · `built` (code + tests) · `validated` (real-data run recorded) · `done` (validated + docs).

| ID | Item | Status | Commits | Validated on |
|---|---|---|---|---|
| W0.1 | Stage status contract, no silent failure | built | develop: `6af3383` | unit (25 tests), stub DAG, Docker fixture run; real-data run pending W0.4 |
| W0.2 | Pinned versions and full provenance | built | develop: `8deea53` | unit (48 tests total), stub DAG local + docker engine (provenance.json valid, 9/9 images resolved), two independent image builds identical (453 packages); real-data run pending W0.4 |
| W0.3 | Test framework and CI | built | develop: `8b25b82`, `a0a7d7d` | 127 unit tests (93 % coverage of bin/ + fungiforge/), 28 nf-test cases, lint 0 warnings, shellcheck clean, stub DAG, Docker fixture run; CI green on develop |
| W0.4 | Cluster development deployment (validate items on real data without touching production) | todo | | |
| W1.1 | Validation suite, controls, benchmarks | built | develop: `ea8c06b`, `6732256` | 142 unit tests incl. minimap2 integration; expected files for CEA10/DF-005; controls chosen (PRJNA1494181) — read-level verification, benchmark runs and control pipeline runs pending on the dev deployment |
| W2.1 | QC and verdict gating in the DAG | built | develop: W2.1 commit | stub DAG + nf-test (forced non_fungal, forced QC fail, --force_all), unit test for the master column; real case (DF-003, bacterial) pending on the dev deployment |
| W2.2 | Read QC: trimmed reads, k-mer profile, in-pipeline triage | todo | | |
| W2.3 | Identification: secondary loci, real concordance, species-aware BUSCO | todo | | |
| W2.4 | Resistance: FungAMR-derived panel, read-level genotyping | todo | | |
| W2.5 | Annotation: GeneMark, species-aware training, InterProScan, eggNOG | todo | | |
| W2.6 | Auxiliary biology: secretome, CAZymes, effectors, ploidy, mating type | todo | | |
| W2.7 | Organelle assembly and annotation | todo | | |
| W2.8 | Mobile elements and mycovirus, made functional | todo | | |
| W2.9 | Novelty with a staged reference genome set | todo | | |
| W3.1 | Phylogenomics and clonality | todo | | |
| W3.2 | Cohort BGC families and mycotoxin flags | todo | | |
| W3.3 | Reporting: per-isolate HTML, cohort MultiQC, merged master table, schema | todo | | |
| W4.1 | Engineering hygiene and documentation pass | todo | | |
| W4.2 | Database fetch, installer and site config updates for new resources | todo | | |
| W5.1 | Cluster upgrade to v0.2.0 and re-validation | todo | | |

## Limitations captured (inventory)

Everything the assessment and the code audit found, mapped to the item that closes it.

| # | Limitation | Evidence | Closed by |
|---|---|---|---|
| L1 | Blanket `\|\| true` on medaka, RepeatModeler/Masker, QUAST, BUSCO, funannotate predict/annotate, barrnap, ITSx, vsearch, sourmash, skani, geNomad, minimap2 | modules stage01–12 | W0.1 |
| L2 | Floating image tags (`latest` ×5), `>=` conda pins, provenance records six tools only — **resolved by W0.2** (versioned tags + digests, explicit conda lock, DB_CHECK, provenance.json with per-tool versions and image sha256) | base.config, env/*.yml, make_provenance.py | W0.2 |
| L3 | No unit tests, no CI; stub DAG tests plumbing only | repo | W0.3 |
| L4 | Resistance detection unproven on a positive; no read-level genotyping, no zygosity; panel has ~35 rows though FungAMR has thousands; species coverage is thin (no *A. flavus* cyp51C, *Cryptococcus*) | af_resistance_panel.tsv, af_resistance.py | W1.1, W2.4 |
| L5 | Non-fungal or QC-failed isolates still run every heavy stage | subworkflow | W2.1 |
| L6 | fastp output discarded; assembly uses untrimmed reads | stage01 | W2.2 |
| L7 | No genome-size, coverage, heterozygosity or ploidy estimate from reads | — | W2.2, W2.6 |
| L8 | Read-level triage is an external script, not a stage | bin/triage_reads.sh | W2.2 |
| L9 | "GCPSR multi-locus concordance" is ITS + ANI only; no TEF1α/RPB2/BenA/CaM, which ITS cannot resolve in *Aspergillus*, *Penicillium*, *Fusarium* | id_classify.py | W2.3 |
| L10 | BUSCO `auto` lineage cannot be species-aware because QC runs before identification | subworkflow order | W2.3 |
| L11 | MLST claimed in stage 08 comments but not implemented | stage08 | W2.3 |
| L12 | GeneMark unlicensed; InterProScan only consumes a pre-computed XML; eggNOG not staged; no RNA-seq evidence path | stage07 | W2.5 |
| L13 | Stage 13 is keyword tallies: secretome, CAZyme, virulence, mating type all effectively empty; ploidy `NA` | extras.py | W2.6 |
| L14 | Mitochondrial split is an empty placeholder file | stage04 | W2.7 |
| L15 | geNomad database is not in the fetch script, so stage 10 never ran; it targets the (empty) mito file; `te_percent` always NA; RVDB-prot is fetched but no stage reads it (found by W0.2's database audit) | stage10, fetch_references.sh | W2.8 |
| L16 | `refseq_fungi_genomes` never staged, so genome-ANI novelty never ran; novelty is ITS-distance only | stage12, fetch_references.sh | W2.9 |
| L17 | No cross-isolate phylogeny, SNP distances or clonality | — | W3.1 |
| L26 | Promoter tandem repeats (TR34/TR46) are called from the assembly only; a duplication absent from the assembly is invisible, and a genotype cannot be confirmed from reads. Found 2026-09-15: the reads deposited as C87 (ERR10820709 + ERR9791656) lack TR34 (0/84 ONT, 1/76 Illumina spanning reads) although the published C87 assembly GCA_949125185 carries it — the detector reports TR34 correctly on that assembly (now a fixture), so the accessions, not the detector, are wrong; a read-verifiable positive control is still needed | cyp51a_TR.py, af_resistance.py | W1.1, W2.4 |
| L27 | Medaka left ~0.3 % residual error at the C87 cyp51A locus (Polypolish changed 0.28–0.31 % of positions per contig): the medaka model is taken from `--dorado_model`, which does not match reads basecalled elsewhere (public data) | stage03 | W2.2 |
| L29 | Stage 05's optional QUAST call could never succeed (QUAST is in none of the images), so every run's assembly-QC status was `partial` and `stages_failed` never read `none`; the call is removed (statistics come from `assembly_qc.py`; the reference benchmark is `bin/benchmark_assembly.py`) | stage05 | W1.1 (done) |
| L28 | `polish.json` records `polypolish_bp_changed` as empty (the Polypolish log is not parsed); C87 produced no ITS at all (`its_present=false`), so identification fell back to sourmash at medium confidence — cause found 2026-09-16: barrnap fails with ENOSPC when the node's /tmp tmpfs is full (left full by crashed annotation tasks), and v0.1.x masked it with `\|\| true`; every task now gets its own work-dir /tmp (v0.1.3 / develop), and identify results of the production batch must be re-checked for silently missing ITS | stage03b, stage08 | W2.2, W2.3 |
| L18 | No BGC families across isolates, no mycotoxin cluster flags | bgc_summary.py | W3.2 |
| L19 | No per-isolate HTML report, no MultiQC, no merged `master_fungi.tsv` produced by the pipeline (Layer 2 gathers rows itself), no output schema | stage14, analysis/ | W3.3 |
| L20 | medaka is given the Dorado model name; a mismatch is hidden by `\|\| true` and falls back to the unpolished assembly silently | stage03 | W0.1, W2.2 |
| L21 | Samplesheet and metadata are not validated (compartment/facility/season free text) | main.nf, cli.py | W4.1 |
| L22 | Manual describes mito, mycovirus and extras as working | docs/manual | W4.1 |
| L23 | Layer-2 R scripts are skeletons that assume columns not yet populated (`has_af_resistance`, `azole_R`, `echino_R`) | analysis/scripts | W3.3, W4.1 |
| L24 | Fetch script downloads through Docker for antiSMASH/funannotate/eggNOG; needs Apptainer path | fetch_references.sh | W4.2 |
| L25 | No chromosome-level scaffolding or telomere check for long-read assemblies; no RNA-seq path | — | out of scope for v0.2 (recorded) |

## Work items

### Phase 0 — Foundations

#### W0.1 Stage status contract
**Goal.** Make silent failure impossible. **Design.** A shared `bin/ff_status.py` helper writes `{status, tools: {name: {version, exit}}, note}` into every stage JSON; modules capture each tool's exit code explicitly; a tool that is optional (e.g. compleasm fallback to BUSCO) records `skipped` with the reason; anything else non-zero is `failed` and the task exits non-zero unless the stage is declared best-effort (BGC, extras), in which case `failed` is recorded and the run continues. `make_report.py` adds `stages_failed` (semicolon list) to the master table. **Validation.** Unit tests for the helper; a stub-run injection test that forces one tool to fail per stage and asserts the JSON status and the master column; CEA10 and DF-005 reruns must show all `ok`.

#### W0.2 Pinned versions and full provenance
**Goal.** Bit-for-bit reproducibility. **Design.** Replace every `latest` with a versioned tag and record the digest in `conf/base.config`; pin conda specs exactly (`==`) with a `conda-lock` file; `make_provenance.py` gathers all tool versions (from each container), image digests, database releases from `MANIFEST.tsv`, git commit, Nextflow version, and the effective parameters; a `--check-db` launch step verifies database completion markers and versions and fails early. **Validation.** Rebuilt images produce identical version tables on two hosts; provenance JSON validated against a schema in the test suite.

**As built (2026-09-15).** Images: `staphb/flye:2.9.6`, `staphb/spades:4.3.0`, `staphb/medaka:2.2.2`, `staphb/polypolish:0.6.1-bwa`, `ezlabgva/busco:v5.7.1_cv1`, `dfam/tetools:2.00` (the validated `latest`; `1.99` ships RepeatModeler 2.0.8/RepeatMasker 4.2.3, so it was rejected after checking), `nextgenusfs/funannotate@sha256:2eadfeb4…` (no versioned tag for the validated 1.8.17.dev197 build); each tag's manifest digest is recorded in `conf/base.config`, and `bin/fetch_references.sh images` pre-pulls exactly these. Conda: `env/base.yml` pinned `==` to the versions of the validated image; `env/base.linux-64.lock` is an explicit lock exported from a fresh build of those pins (the Mac's 0.1.0 image turned out to lack kraken2 — the lock check caught it — so the lock comes from a rebuild that includes it); `env/Dockerfile` and `env/fungiforge.def` install from the lock (`--build-arg ENV_SPEC=env/base.yml` bootstraps a new lock); `bin/lock_env.sh` regenerates and cross-checks. Provenance: `ff_version` in the status helper records each tool's version into the stage JSON (`tools[].version`, `versions`); `DB_CHECK` (stage 00, run-level, `bin/check_databases.py`) gates every isolate's first stage and fails the run with one list of missing databases unless `--allow_missing_db`; `PROVENANCE` (stage 15, run-level, `bin/make_provenance.py aggregate`) writes `pipeline_info/provenance.json` and validates it against `fungiforge/resources/provenance.schema.json`; `workflow.onComplete` adds image identities (sha256 of `.sif`/cached `.img`, or docker repo digests; `--provenance_hash_images`, `--image_cache_dir`). Basecall (stage 00) was also brought onto the status contract (it had been missed by W0.1). Validation done: 48 unit tests (`ff_version` end to end incl. missing/failing commands; DB checker; aggregate/images/validate; schema semantics); stub DAG 46/46 with `local` and with the `docker` engine (all nine images resolved with repo digests; document schema-valid); every version command exercised in the rebuilt image; two independent builds (from `base.yml` and from the lock) have identical 453-package sets. Still pending: real-data run on the dev deployment (W0.4), and the apptainer cache path on the cluster. Note for the record: tool-reported strings can differ from package versions (chopper 0.9.0 prints `0.8.0`, barrnap 1.10.6 prints `1.10.5`); the lock is authoritative.

#### W0.4 Cluster development deployment
**Goal.** Validate each item on real data while production keeps running. **Design.** A second, independent install at `/hpc/opt/fungiforge-dev` (`hpc_install.sh --prefix /hpc/opt/fungiforge-dev --ref develop --db /hpc/data/fungiforge --group <group>` with its own images and site config, no profile.d hook), run directories under `~/runs/dev/`. Items are validated there on CEA10 and DF-005 (and C87 once available) before being marked `validated`. **Validation.** The dev install runs the stub DAG and CEA10 end to end with the develop branch.

#### W0.3 Test framework and CI
**Goal.** Every later item ships with tests that run automatically. **Design.** `pytest` for all `bin/*.py` with fixtures under `test/fixtures/`; nf-test cases for the routing logic and each module's stub; GitHub Actions running unit tests, the stub DAG on the three-mode fixture, `nextflow lint`, and shellcheck; a nightly optional job running the CEA10 subsample through Docker. **Validation.** CI green on `develop`; coverage report for `bin/`.

**As built (2026-09-15).** Unit: 127 tests over every `bin/*.py` (13 new test files; synthetic GenBank/protein/promoter/Kraken2 fixtures, `test/fixtures/kraken2/`, shared `conftest.py`) and the CLI; scripts are driven as subprocesses with subprocess coverage (`.coveragerc`, `sitecustomize`), 93 % of `bin/` + `fungiforge/`, CI fails below 85 %. The tests found and fixed three real defects: `id_classify` never dropped `Genus_sp` hits (the `_sp$` check ran after underscores were replaced, so "Aspergillus sp" was reported as a species-level call), a panel row whose organism field contained "spp." applied to every genus in `af_resistance`, and a duplicated accession in one FungAMR cell double-counted its mutations; `make_report` no longer crashes on a `top_species` entry without `pct`. Static: `nextflow lint` at zero warnings — `Channel` → `channel`, explicit closure parameters, unused-parameter prefixes, and every task now calls the helpers from `bin/` on PATH (`source ff_status.sh`, `make_report.py …`) instead of `${projectDir}/bin/…` (verified inside the real containers by the Docker fixture run; all `bin/*.py` are executable); shellcheck clean at warning level on every shell script (directives only where the dynamic `source` is intended). nf-test 0.9.5 (`nf-test.config`, `test/nf-test/`): six routing cases on `parse_row`/`resolve_path`, one stub test per module (19), and three whole-DAG cases (46 tasks; skip flags remove 12; `--help` launches nothing) — 28 cases. First CI run failed on the two Nextflow jobs because the per-label memory requests (16 GB) exceed a hosted runner; `process.resourceLimits` in `conf/test.config` now clamps the machinery check, and run 34944666364 is green on all five jobs. `test/run_docker_fixture.sh` runs the fixture through the real images: DB_CHECK and READ_QC complete in-container with versions recorded, and Flye/SPAdes fail on the synthetic reads with the failure recorded by the contract (the expected outcome; it is the in-container check the hosted CI cannot do). CI (`.github/workflows/ci.yml`): unit (Python 3.11/3.12, coverage artifact), lint, stub DAG (validates `provenance.json`), nf-test. **Deviation from the design:** the nightly CEA10 Docker job is not feasible on a hosted runner (~40 GB databases, ~30 GB images); real-data validation stays on the development deployment (W0.4/W1.1) and the Docker fixture run is documented as a local check.

### Phase 1 — Validation

#### W1.1 Validation suite, controls and benchmarks
**Goal.** Demonstrated, not designed. **Design.** Controls: C87 (TR34/L98H), a *C. auris* FKS1 S639F isolate and a *C. albicans* ERG11 isolate from the SRA (accessions recorded here when chosen), the wild-type siblings C6 and E142; reference benchmarks: CEA10 against its published assembly and DF-005 against *A. flavus* NRRL 3357 using QUAST (misassemblies, NGA50, indels per 100 kb) and BUSCO; `test/expected/<sample>.json` with tolerances; `bin/validate_run.py` producing a pass/fail table. **Validation.** All controls call their known genotype at high confidence; benchmarks within tolerance; results recorded in `docs/validation.md`.

**As built (2026-09-16).** `bin/validate_run.py` (also `fungiforge validate`) compares a run with `test/expected/<sample>.json` — rules `equals`, `one_of`, `regex`, `min`, `max`, `expected` ± `abs_tol`/`rel_tol`, values read from stage JSONs, the master row or provenance — and prints a pass/fail table (exit 1 on failure); expected files for CEA10 (hybrid) and DF-005 (Illumina-only) written from the recorded v0.1.0 values with explicit tolerances. `bin/benchmark_assembly.py` gives the QUAST-style numbers (aligned fractions, NGA50, breakpoints, mismatches and indels per 100 kb, unaligned contigs) from minimap2 + the cs tag, because QUAST is in none of the images; `fetch_references.sh benchmarks` stages A1163 and NRRL 3357. Controls: the C87 accessions were shown to be a different isolate from the published assembly (see docs/validation.md §3.1), so controls are now taken from PRJNA1494181 (ONT, Dorado hac): two TR34/L98H isolates and a matched wild-type sharing the same cyp51A background; `bin/control_genotype.py` verifies a control's genotype from its reads at a wild-type locus (`test/controls/afum_cyp51A_wt_locus.fa`, derived from the published C87 assembly with one TR unit removed and L98H reverted — identical to AF338659 over its 352-bp promoter and exon structure) with insertion-window and point-mutation checks, and `bin/fetch_controls.sh` downloads, verifies and writes a samplesheet row per control, refusing any whose reads fail the check. Tests: 15 new (rules, table/exit status, SAM walking, call thresholds, spec-locus consistency, minimap2 integration on simulated ONT reads for TR34/L98H and wild-type, benchmark metrics on hand-built and synthetic alignments); minimap2 added to the CI unit job. Also fixed: the phantom QUAST step (L29). *Candida* controls deferred to W2.4. **Pending for `validated`:** the dev-deployment CEA10/DF-005 runs compared with the expected files, the two reference benchmarks, and the three poultry-farm controls fetched, read-verified and run.

### Phase 2 — Correctness of existing stages

#### W2.1 QC and verdict gating in the DAG
**Design.** After stage 04 and 05, isolates with `verdict in (non_fungal, human)` or `qc_pass == false` are routed to a `SKIPPED_REPORT` path that emits the master row with `stages_failed` and the reason, skipping stages 06–13; `--force_all` overrides. **Validation.** Stub test with a forced non-fungal verdict; DF-003 (bacterial) as the real case: row present, no heavy stage executed.

**As built (2026-09-16).** The subworkflow reads each isolate's decontam and assembly-QC JSONs after Stage 05 and stops isolates with verdict `non_fungal`/`human` or `qc_pass == false` (`mixed` and `not_run` pass); a small `GATE` process (`modules/stage05b_gate.nf`) writes `<id>.gate.json` (status `skipped`, reason, the eight stages not run) so REPORT still produces the master row — new append-only column `gate` (`pass` | `skipped(<reason>)`) and `gate:skipped` in `stages_failed`; stages 06–13 receive only isolates that passed (`nuclear_ok` channel), everything downstream inherits the filter through its joins. `--force_all` disables the gate. Test hooks for stub runs: `--stub_verdicts 'ID=non_fungal'` and `--stub_qc_fail 'ID'` make the decontam/QC stubs emit the forcing values; nf-test cases assert the task count (46 − 2×8 + 2), that no fungal stage ran for the stopped isolates, the gate JSONs and the provenance entries, and that `--force_all` restores 46 tasks. Folded in: Medaka now resolves its model from the reads' `basecall_model_version_id` header (Dorado) and otherwise uses `--medaka_model` (default: the `--dorado_model` string) after checking it against `medaka tools list_models` — an unknown model fails the stage instead of polishing with the wrong one (L27); the Polypolish change count is summed from its per-contig log lines into `polish.json` (L28, first half). **Pending for `validated`:** DF-003 (bacterial) through the dev deployment: row present with `gate = skipped(verdict:non_fungal)`, no annotation task; a fungal isolate unchanged.

#### W2.2 Read QC upgrade
**Design.** Route fastp-trimmed reads to assembly and polishing; add k-mer profiling (KMC3 + GenomeScope2) giving genome size, heterozygosity, ploidy hint and expected coverage into `readqc.json`; add an optional in-pipeline read triage (Kraken2 on a subsample) as stage 00b writing the verdict early so W2.1 can gate before assembly; correct medaka model mapping (derive from the Dorado model string, validate against `medaka tools list_models`, fail if unknown). **Validation.** Unit tests for the model mapping and GenomeScope parsing; DF-005 genome size within 10 % of assembly length; CEA10 medaka runs with the correct model (log inspected).

#### W2.3 Identification upgrade
**Design.** Extract secondary loci from the assembly with curated reference sets and HMMs: LSU (D1/D2), TEF1α, RPB2, β-tubulin (BenA), calmodulin (CaM); search each against a curated marker database (UNITE for ITS; per-locus references from RefSeq/GenBank type material, staged by W4.2); true concordance: species call requires ITS and at least one secondary locus to agree, else genus-level with a flag; run BUSCO twice: `fungi_odb10` for the early gate, then the class/order lineage chosen from the species call for the reported completeness; implement MLST via `mlst` for schemes that exist (*A. fumigatus*, *C. albicans*, *C. glabrata*, *C. tropicalis*, *C. krusei*). **Validation.** DF-005 resolves *A. flavus* vs *A. oryzae* on CaM/BenA; CEA10 stays *A. fumigatus*; unit tests on marker extraction.

#### W2.4 Resistance upgrade
**Design.** Generate the panel from FungAMR programmatically (`build_fungamr_refs.py` extended) with species, gene, position, wild-type, mutant, drug class, evidence tier; keep the assembly-level path; add read-level genotyping: map reads (minimap2 for ONT, bwa-mem2 for Illumina) to the species reference gene set, call variants with bcftools, translate to protein changes, report allele frequency and zygosity; confidence combines polish mode, read support and evidence tier; species-specific extras: *A. fumigatus* cyp51A TR (kept), *A. flavus* cyp51C hotspots, *Candida* ERG11/FKS1/FKS2 hotspot regions, *Cryptococcus* FCY1/FUR1/ERG11; copy-number signal for ERG11/CDR1 from coverage. **Validation.** Controls in W1.1 call their genotypes; CEA10/DF-005 remain wild-type; unit tests on panel generation and variant translation.

#### W2.5 Annotation upgrade
**Design.** GeneMark-ES via `--genemark_key` (key installation documented, Appendix A); species-aware `--busco_db` and Augustus species from the stage-08 call; InterProScan as its own containerised stage (07b) with the InterPro data staged by W4.2, feeding `funannotate annotate --iprscan`; eggNOG-mapper database staged; BRAKER3 offered as `--annotator braker3` when RNA-seq is provided (future). **Validation.** DF-005 and CEA10 gene counts and BUSCO-protein completeness compared before and after; Pfam/GO coverage reported.

#### W2.6 Auxiliary biology
**Design.** SignalP 6 (academic licence, Appendix A) and EffectorP 3 for secretome and effectors; dbCAN3 (`run_dbcan`) for CAZymes; nQuire and Smudgeplot on the reads for ploidy and heterozygosity; MAT locus typing by HMM (Pfam MATalpha_HMGbox PF04769, HMG box PF00505 in the MAT1-2-1 context) plus synteny check. **Validation.** *A. fumigatus* CEA10 is MAT1-1 (published); DF-005 haploid; counts non-zero and within published ranges.

#### W2.7 Organelle assembly
**Design.** GetOrganelle (`fungus_mt`) from reads, oatk as alternative for ONT; MITOS2 or MFannot annotation; mito size, gene set, introns, and homoplasmy; feeds stage 10. **Validation.** CEA10 mitogenome ~32 kb with the expected gene set.

#### W2.8 Mobile elements and mycovirus
**Design.** Stage geNomad's database; run on the whole assembly (plasmids/proviruses, EVEs) and the mitogenome; RepeatMasker `.tbl` parsed into `te_percent` and class breakdown; mycovirus screen by DIAMOND blastx of unplaced/low-coverage contigs against RVDB-prot, with the DNA-only caveat stated. **Validation.** `te_percent` populated for both isolates; geNomad status `ok`.

#### W2.9 Novelty
**Design.** Stage a representative fungal genome set (RefSeq representative genomes plus type-strain genomes for genera in the study, via NCBI datasets, listed in the manifest); skani ANI against it; novelty from ANI (< 95 % to nearest = candidate novel species, 95–97 % flagged) combined with ITS distance and the secondary loci. **Validation.** CEA10/DF-005 known; a synthetic divergent genome flagged novel in a unit test.

### Phase 3 — Cohort features

#### W3.1 Phylogenomics and clonality
**Design.** Post-isolate stage: BUSCO single-copy orthologs across isolates → MAFFT → trimAl → IQ-TREE species tree; within same-species clusters, reference-based SNP calling (minimap2/bwa + bcftools, MycoSNP-style filters) → pairwise SNP distances and clonal groups at a documented threshold; outputs to `results/cohort/`. **Validation.** CEA10 vs C87 vs C6/E142 tree topology matches published relationships; unit tests on distance matrices.

#### W3.2 Cohort BGC families and mycotoxin flags
**Design.** BiG-SCAPE over all region GenBanks with MIBiG; families table; mycotoxin flags from KnownClusterBlast hits to curated MIBiG entries (aflatoxin, ochratoxin, fumonisin, gliotoxin, patulin, trichothecene, ergot). **Validation.** DF-005 flags the aflatoxin cluster; CEA10 flags gliotoxin.

#### W3.3 Reporting
**Design.** Per-isolate HTML report (Jinja2) with QC, identification, resistance, BGC and mobile summaries; cohort MultiQC; the pipeline writes a merged `results/master_fungi.tsv` via `collectFile`; a JSON schema for the master table checked in the report stage; Layer-2 scripts updated for new columns. **Validation.** Schema validation in tests; Layer-2 runs end-to-end on the batch output.

### Phase 4 — Engineering and documentation

#### W4.1 Hygiene and documentation pass
Samplesheet and metadata schema validation with clear errors; `nextflow lint` clean; shellcheck clean; manual rewritten to describe actual behaviour with a roadmap section; CHANGELOG; CITATION updated; v0.2.0 release notes.

#### W4.2 Databases, installer, site config
Fetch script gains geNomad, reference genome set, marker-locus references, InterPro data, eggNOG, SignalP model placement; Apptainer path for the container-driven downloads; installer and site config updated for new images and resources; deployment guide updated.

### Phase 5 — Deployment

#### W5.1 Cluster upgrade and re-validation
After production runs complete: merge `develop` → `main`, tag `v0.2.0`, `hpc_install.sh --ref v0.2.0 --rebuild-images`, rerun CEA10, DF-005 and all controls, run `validate_run.py`, update the deployment document and `docs/validation.md`.

## Appendix A — Licensed software

**GeneMark-ES/ET/EP+** (Georgia Tech). Free for academic use. Request at `http://exon.gatech.edu/GeneMark/license_download.cgi`: choose "GeneMark-ES/ET/EP+" for LINUX 64, fill in name, institution and email, accept the licence; download `gm_key_64.gz`, then `gunzip gm_key_64.gz && cp gm_key_64 ~/.gm_key` on the machine that runs annotation, or point `--genemark_key` at it. Keys expire roughly yearly; renew the same way. The GeneMark binaries themselves are bundled in the funannotate image; only the key is needed.

**SignalP 6.0** (DTU Health Tech). Free academic licence: `https://services.healthtech.dtu.dk/services/SignalP-6.0/`, "Downloads", register with an academic email, download the `signalp-6.0h.fast.tar.gz` package, and place it where W2.6's build expects it (documented there). EffectorP 3, dbCAN3, nQuire, Smudgeplot, GetOrganelle, MITOS2, geNomad, BiG-SCAPE and InterProScan are free and need no registration.

## Change log of this document

- 2026-09-15 — created; limitation inventory L1–L25; items W0.1–W5.1 defined; all `todo`.
- 2026-09-16 — W2.1 built: the gate after assembly QC (GATE process, `gate` master column, `--force_all`, stub hooks, nf-test cases); Medaka model resolved from read headers with validation (L27); Polypolish change count recorded (L28 first half).
- 2026-09-16 — W1.1 built: validate_run / benchmark_assembly / control_genotype tools with tests, expected files for CEA10 and DF-005, PRJNA1494181 controls chosen with read-level verification, fetch_controls.sh, docs/validation.md; phantom QUAST step removed (L29); *Candida* controls deferred to W2.4.
- 2026-09-16 — Node /tmp tmpfs incident, third round: v0.1.1's TMPDIR export did not reach funannotate-p2g.py (hard-coded /tmp); v0.1.2 bound a work-dir directory over /tmp for the annotate task; the full tmpfs then broke barrnap (ENOSPC, exit 28) in the dev CEA10 run — visible only because of the W0.1 contract. v0.1.3 / develop: EVERY containerised task gets its own /tmp (process-level `beforeScript 'mkdir -p tmp'` + `containerOptions` bind) under docker/singularity/apptainer; verified with Docker on plain and antiSMASH (entrypoint-cleared) tasks.
- 2026-09-15 — C87 positive control run (production v0.1.0, 15/15 stages): cyp51A Y121F + T289A + G448S at high confidence, promoter wild-type. Read-level check showed the deposited reads genuinely lack TR34 while the published assembly has it; the detector called TR34 on the published assembly (fixture `test/fixtures/cyp51A/`). Added L26–L28. The production batch also hit a node-tmpfs overflow in annotation (v0.1.1 hotfix on main, ported here).
- 2026-09-15 — W0.3 built: 127 unit tests over every helper and the CLI (subprocess coverage 93 %), `nextflow lint` zero warnings (tasks call `bin/` helpers via PATH, no `projectDir` in processes), shellcheck clean, nf-test suite (routing + 19 module stubs + DAG, 28 cases), Docker fixture script, GitHub Actions CI. Real bugs fixed on the way: `id_classify` species-level calls on `Genus_sp` hits, `af_resistance` "spp." matching every genus, duplicate FungAMR accessions double-counted. Nightly Docker CEA10 job dropped as infeasible on hosted runners (documented).
- 2026-09-15 — W0.2 built: pinned image tags + digests, exact conda pins + explicit lock (image builds install from it), `ff_version`, `DB_CHECK` and `PROVENANCE` run-level stages, schema-validated `provenance.json` with image sha256/digests, basecall stage onto the status contract, 23 new unit tests, stub DAG with and without a container engine, lock reproducibility verified by two builds. Found on the way: the local 0.1.0 image lacked kraken2 (lock check), `dfam/tetools:1.99` is older than the validated `latest` (now `2.00`), RVDB is fetched but unused (noted under L15).
- 2026-09-15 — W0.1 built: `bin/ff_status.sh` + `bin/ff_status.py`, all 15 modules converted (no `|| true` masking a tool), assemble/medaka status JSONs, `stages_failed` master column, 25 unit tests, Docker fixture run exercised ok and failed paths. Found and fixed in the process: a bare optional `ff_run` tripped `bash -ue` errexit (now returns 0 and sets `$FF_RC`). Added W0.4 (cluster dev deployment) because real-data validation must not touch the production install.
