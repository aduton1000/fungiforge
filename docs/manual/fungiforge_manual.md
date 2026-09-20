```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

# 1. Introduction

**FungiForge** is a **Nextflow DSL2** pipeline that turns sequencing reads from a single
fungal isolate — **Oxford Nanopore (ONT) long reads, Illumina short reads, or both (hybrid)** —
into a **de-novo assembly, a eukaryotic gene annotation, a species
identification, an antifungal-resistance genotype, a mobile-element/mycovirus inventory, a
biosynthetic-gene-cluster (BGC) catalogue, and a novelty verdict**, then merges everything
into one per-isolate report and a single row of a fixed-schema **master table**. A second,
independent **R analysis layer** consumes those master rows across many isolates to build the
air / human / surface **One Health** comparative story.

FungiForge is the **fungal sibling of the *forge* family** (CaptureForge designs capture
panels, CallForge analyses targeted-capture germline callsets, MethylForge analyses
methylomes). It follows the same engineering stance: **stand on validated upstream tools,
add only the bespoke layers the science genuinely needs, and change inputs and configuration
— never code — to move to a new organism, facility, or study.**

**Its immediate use** is the fungal arm of an airborne-AMR One Health study: culturing and
sequencing fungi recovered from **air, human, and surface** compartments across abattoirs,
clinics, markets and community settings, to ask which settings concentrate azole-resistant
*Aspergillus fumigatus*, echinocandin-resistant *Candida/Nakaseomyces*, priority pathogens
such as *Candida auris*, and candidate novel species or novel natural-product clusters. But
the pipeline is **organism-agnostic**: any cultured fungal isolate with ONT and/or Illumina
reads runs through it.

**Who it is for.** Genomics-core and bioinformatics staff running fungal isolate genomics
who need a reproducible, container-pinned workflow that is honest about the two structural
realities of fungal genomics: **(1)** ONT homopolymer indels fall exactly where
antifungal-resistance point mutations live, so resistance calls on an ONT-only assembly are
**provisional until hybrid-polished** (hybrid and Illumina-only assemblies are high-confidence,
having no homopolymer-indel problem); and **(2)** there is **no ResFinder for fungi and no
GTDB for fungi** — resistance and identification must be done with bespoke, curated logic
rather than a single push-button database.

> **Status note used throughout this manual.** FungiForge is implemented and validated
> **end-to-end on the synthetic `test` fixture** (all 15 stages wire and run under
> `-profile test`). Stage implementations are being exercised on a first real *A. fumigatus*
> isolate. Where a statement is a property of the *machinery* rather than a real biological
> result, the text says so. This manual documents **what the code actually does** and flags
> what is a placeholder or milestone-pending honestly (see §12–§14).

# 2. Concepts

This section covers the fungal-genomics background a FungiForge user needs to read the
outputs correctly.

**Three input modes — and why polishing matters for resistance.**
FungiForge assembles from whatever reads an isolate has, routed automatically by
`meta.assembly_mode` (set from the samplesheet):

- **ONT-only (`longread`)** — Nanopore long reads span repeats and resolve a contiguous fungal
  genome (typically 30–50 Mb) that short reads alone fragment (Flye). Their weakness is
  **systematic error in homopolymer runs** — a run of six `A`s may be read as five or seven —
  producing **false indels and frameshifts**. Antifungal-resistance calling reads single-residue
  changes in genes like *cyp51A*, *ERG11* and *FKS1*, so an uncorrected homopolymer error can
  *invent* a resistance mutation or *erase* a real one. **Stage 03a Medaka** (an ONT neural
  consensus) is always applied, but calls remain `provisional_ont_only`.
- **Hybrid (`longread` + Illumina)** — after Medaka, **Stage 03b Polypolish** uses the short reads
  to correct the residual indels Medaka cannot. This is the highest-quality assembly and yields
  `high`-confidence resistance calls.
- **Illumina-only (`shortread`)** — no ONT at all: **Stage 02b SPAdes** assembles directly from the
  short reads. There is no homopolymer-indel problem, so Medaka/Polypolish are skipped and calls
  are `high`-confidence. (The trade-off is contiguity: short-read-only assemblies are more
  fragmented than long-read ones, which can affect repeat-spanning features such as the *cyp51A*
  promoter TR locus.)

Every isolate carries a **polish mode** (`hybrid`, `illumina_only`, `ont_only`, or
`hybrid_failed_ont_fallback`) and a **resistance-confidence flag**; ONT-only calls are labelled
`provisional_ont_only` and never presented as clean.

**ITS/UNITE barcoding and GCPSR multi-locus identification — no GTDB for fungi.** Bacteria
have GTDB-Tk and a genome-taxonomy backbone; fungi do not. The universal fungal barcode is
the **ITS region** (ITS1–5.8S–ITS2, between the SSU and LSU rRNA genes), classified against
the **UNITE** reference database of curated species hypotheses (SH). ITS alone is often only
genus-resolving, so FungiForge combines it with **genome-level average nucleotide identity
(ANI)** (sourmash/skani vs type-strain sketches) and **MLST**, and reports a call under the
principle of **GCPSR — Genealogical Concordance Phylogenetic Species Recognition**: a species
is asserted with high confidence only when independent loci (ITS *and* genome ANI) agree. The
result is always a species **plus a confidence and a method**, never a bare label. Guidance
thresholds: ITS identity ≥98.5% → species-level; 94–98.5% → genus/uncertain; <94% →
candidate novel.

**Why antifungal resistance is bespoke: a curated panel + FungAMR references + a structural
TR detector.** No curated, sequence-bearing "ResFinder for fungi" exists. FungiForge builds
resistance calling from three parts:

1. **A curated panel** (`af_resistance_panel.tsv`) — an interpretation overlay listing, per
   gene and per species, the drug class, the *mechanism* (substitution / loss-of-function /
   gain-of-function / over-expression / promoter tandem-repeat), the hotspot residue numbers,
   and the known resistance mutations, each with a literature source.
2. **FungAMR reference sequences** — FungAMR (Bédard/Landry *et al.*, *Nat Microbiol* 2025)
   ships a mutation *table* with no sequences; `bin/build_fungamr_refs.py` fetches the
   UniProt protein for every `ref_seq_uniprot_accession` and writes a species-labelled
   reference FASTA that the caller aligns the isolate's proteins against, so hotspot residue
   numbering is anchored to the correct reference frame.
3. **A dedicated structural detector for the *A. fumigatus* cyp51A promoter tandem repeats.**

**cyp51A TR34/TR46 promoter structure.** The dominant *environmental* azole-resistance
mechanism in *A. fumigatus* is **not** a point mutation but a **tandem repeat (TR) inserted
into the cyp51A promoter** — a duplication of a ~34 bp (TR34) or ~46 bp (TR46) unit in the
region between the two SrbA transcription-factor binding sites, which up-regulates the azole
target. TR34 canonically pairs with the coding mutation **L98H**, and TR46 with **Y121F +
T289A**. Because a TR is a structural insertion, it is invisible to any substitution caller;
`bin/cyp51a_TR.py` locates cyp51A in the annotation, extracts the ~600 bp upstream promoter
(strand-aware), and scans for a tandem array, reporting the unit length, copy number and TR
type (TR34 / TR46 / TR46³ / TR53).

**Biosynthetic gene clusters (BGCs).** Fungi encode clustered secondary-metabolite pathways —
non-ribosomal peptide synthetases (NRPS), polyketide synthases (PKS), terpene synthases,
RiPPs, and hybrids — that produce mycotoxins, antibiotics and virulence factors.
**fungiSMASH (antiSMASH 8, fungal mode)** detects them per isolate; whether a cluster is a
known MIBiG compound or novel is resolved across all isolates by BiG-SCAPE in the Layer-2
comparative step.

**Novelty by ANI/ITS distance.** A candidate novel species is flagged when **genome ANI to
the closest reference genome is <95%** (the strong signal, when a reference genome set is
staged) **or**, as a fast fallback, when the **best ITS identity is <98.5%**. Both verdicts
are reported with the GCPSR/polyphasic caveat — a genome distance is a hypothesis, confirmed
only by multi-locus phylogeny and a formal polyphasic description.

# 3. Pipeline architecture

FungiForge Layer 1 runs **15 Nextflow stages** (`modules/stage00…14`), chained by
`subworkflows/fungiforge.nf` and launched from `main.nf`. Each stage runs in its own pinned
container (or conda env), publishes its outputs to `results/<sample>/<NN>_<stage>/`, and
emits a small **`result.json`**; Stage 14 merges every stage's JSON into a per-isolate HTML
report and one row of the master table. The diagram below is the authoritative map: solid
edges are the assembly spine, dotted edges feed a downstream stage a specific artifact,
dashed nodes are optional/skippable, and the amber cylinders are the reference databases
staged under `--data_dir`.

![FungiForge Layer-1 pipeline DAG — the 15 stages colour-coded by group; the optional Dorado basecall and the four skippable feature-calling stages are dashed; external reference databases (amber) feed the ID/resistance/BGC/novelty stages; the antifungal-resistance stage (red, ★) is the bespoke centerpiece; everything converges on Stage 14, which emits `master_fungi.tsv` — the handoff to the Layer-2 One Health analysis.](figures/architecture.png){width=100%}

**Stage walkthrough (purpose · key tool/container · inputs → outputs).**

| # | Stage | Key tool(s) / container | Inputs → key output(s) |
|:--|:----------|:--------------------|:-----------------------------|
| 00 | Basecall *(optional)* | Dorado (native host binary / base image) | pod5 dir → `*.ont.fastq.gz` |
| 01 | Read QC & filter | NanoPlot, chopper, fastp / base image | raw reads → `*.ont.filt.fastq.gz`, `readqc.json` (mode-aware) |
| 02 | Assembly (long-read) | Flye \| Canu \| Raven / `staphb/flye` | filtered ONT → `*.assembly.fasta` |
| 02b | Assembly (short-read) | SPAdes \| MEGAHIT / `staphb/spades` | Illumina-only → `*.assembly.fasta` |
| 02c | Haplotig purging | purge_dups + minimap2 / base image | long-read draft → `*.purged.fasta`, `purge.json` |
| 03a | ONT polish | Medaka / `staphb/medaka` | draft + ONT → `*.medaka.fasta` (long-read isolates) |
| 03b | Short-read polish (hybrid) / Illumina-only passthrough | bwa + Polypolish / `staphb/polypolish` | Medaka + Illumina → `*.polished.fasta`, `polish.json`; Illumina-only passes through (`mode=illumina_only`) |
| 04 | Decontam + organelle split | Kraken2 per contig; mitochondrial contigs separated by core-gene tblastn / base image | polished → `*.nuclear.fasta`, `*.mito.fasta`, `decontam.json` |
| 04b | Mitochondrial genome | core genes, rRNAs, introns, circularity, copy number, heteroplasmy / base image | mito + reads → `organelle.json`, `*.mito.gff` |
| 05 | Assembly QC + completeness | compleasm/BUSCO + contiguity / `ezlabgva/busco` | nuclear → `assemblyqc.json`, `*.busco_sc.faa` (single-copy proteins for the cohort tree) |
| 05b | Gate | records why an isolate is stopped (non-fungal verdict, failed QC) / base image | → `gate.json` |
| 06 | Repeat model + soft-mask | RepeatModeler2, RepeatMasker / `dfam/tetools` | nuclear → `*.masked.fasta`, `*.telib.fasta`, `repeat.json` |
| 07a | Gene prediction | Funannotate predict (GeneMark-ES when licensed; species-aware training) / `nextgenusfs/funannotate` | masked + species → `predict_results/`, `predict.json` |
| 07b | eggNOG-mapper | emapper vs eggNOG 5 / `eggnog-mapper` image | proteins → `*.emapper.annotations` |
| 07c | InterProScan *(opt-in)* | InterProScan 5 / `interpro/interproscan` | proteins → `*.iprscan.xml` |
| 07 | Functional annotation | Funannotate annotate / `nextgenusfs/funannotate` | prediction + eggNOG + InterPro → `*.proteins.faa`, `*.gbk`, `annotate.json` |
| 08 | Identification (multi-locus concordance) | ITSx, barrnap, vsearch/UNITE, tblastn/blastn vs type-material sets, mlst, sourmash / base image | nuclear → `*.species.txt`, `*.markers.fasta`, `*.busco_lineage.txt`, `identify.json` |
| 08b | Species-aware BUSCO | compleasm/BUSCO with the lineage from stage 08 / BUSCO image | nuclear + lineage → `busco_lineage.json` |
| 09 | **Antifungal resistance ★** | `af_resistance.py` + `cyp51a_TR.py` / base image | proteins + species + nuclear + GBK → `resistance.json` |
| 10 | Mobile & repeat elements | RepeatMasker table, geNomad, RVDB screen, mito homing endonucleases / base image | TE library + assembly + mito → `mobile.json` |
| 11 | Biosynthetic gene clusters | fungiSMASH (antiSMASH 8) / `antismash/standalone:8.0.0` | GBK → `bgc.json`, `*.regions.gbk` |
| 12 | Novelty | skani vs the reference genome set + ITS/secondary loci / base image | nuclear + ID → `novelty.json` |
| 13 | Eukaryote extras | SignalP 6 (site image), EffectorP 3, dbCAN, PHI-base, Pfam MAT locus, nQuire / base image | proteins + GBK + read BAM → `extras.json` |
| 14 | Report + master row | `make_report.py` (Jinja2 template) / base image | all `result.json` → `*.report.html`, `*.master.tsv` |
| 16 | Cohort phylogenomics *(run-level)* | skani clusters, BUSCO supermatrix + IQ-TREE 3, core SNP distances, BiG-SCAPE families / base image | every passing isolate → `cohort/` |
| 17 | Cohort summary *(run-level)* | merged + validated master table, cohort report, MultiQC content / base image | all master rows → `04_summary/` |

**Nextflow / profile model.** `main.nf` parses the sample sheet and calls the single
subworkflow `FUNGIFORGE`, which wires the 15 modules. Behaviour is set by **two composable
profile groups** on the command line — an **execution** profile (`local` or `hpc_slurm`) plus
a **packaging** profile (`conda`, `mamba`, `docker`, `singularity`, or `apptainer`).
`conf/base.config` maps each stage's resource **label** to its container image / conda env and
its CPU/memory/time request; Nextflow picks the container when a container profile is active,
otherwise the conda env. Heavy stages use their upstream community images; the light stages
share the `aduton1000/fungiforge:0.1.0` base image.

**The two-layer design and the master-table contract.** Layer 1 is per-isolate genomics;
**Layer 2 (`analysis/`)** is a set of R objective scripts that read only the master table.
The handoff is a **fixed-schema TSV** so the two layers stay decoupled — Layer 2 can be
re-run independently, and the schema deliberately mirrors the bacterial study's
`master_results` table so the fungal and bacterial One Health stories are directly
comparable.

![Two-layer architecture and the master-table handoff contract. Layer 1 (Nextflow) emits one `master_fungi.tsv` row per isolate to `results/04_summary/`; Layer 2 (R) merges those rows and runs eight One Health objectives.](figures/two_layer.png){width=92%}

# 4. Installation

FungiForge needs **Nextflow (≥23.10)** with **Java 17+**, plus either a container engine
(Docker on a workstation, Apptainer/Singularity on HPC) or conda/mamba. You do **not** install
the bioinformatics tools by hand — each stage declares its own pinned container or `env/*.yml`,
and Nextflow provisions it on first use.

## 4.1 Local install (workstation / laptop)

```bash
# 1. Nextflow (Java 17+ required)
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
#    or: mamba create -n nextflow -c bioconda 'nextflow>=23.10'

# 2. obtain FungiForge and enter the repository
cd /path/to/fungiforge

# 3. install the small Python front-end + bespoke callers (stdlib-light)
pip install -e .          # registers the `fungiforge` command; also usable in a conda env
```

**macOS (Apple Silicon) — the amd64 emulation caveat.** All heavy stage images are built
**`linux/amd64`**; there are no `osx-arm64` builds for tools like Funannotate or antiSMASH.
On Apple Silicon the `docker` profile therefore runs them under **Rosetta/QEMU emulation** —
correct but slow. The pipeline is configured for this: `docker.runOptions` sets
`--platform linux/amd64`, and heavy/slow options default **off** for local runs
(`--run_interproscan false`, `--genome_id false`). Two consequences to plan around:

- Annotation (Stage 07) and BGC detection (Stage 11) can take **many hours** per isolate under
  emulation — budget accordingly, or run them on the HPC.
- **Dorado (Stage 00) runs as a native arm64 host binary, not in a container**, so basecalling
  stays fast on Apple Silicon. The light conda envs that have no `osx-arm64` build can be created with
  `CONDA_SUBDIR=osx-64 mamba env create -f env/<name>.yml` (noted in `env/readqc.yml`).

**The container image.** The light stages use `aduton1000/fungiforge:0.1.0`; build it once
(native, no emulation) from the repo:

```bash
docker build -f env/Dockerfile -t aduton1000/fungiforge:0.1.0 .
```

## 4.2 Shared HPC install (SLURM + Apptainer/Singularity)

The production route is the cluster (native linux-64, no emulation). `docs/hpc_deployment.md`
is the authoritative procedure; in brief:

```bash
# 1. code + tools (nextflow >=23.10, java 17+, apptainer or singularity)
git clone https://github.com/aduton1000/fungiforge && cd fungiforge

# 2. build the base image once (heavy public images pull automatically at run time)
apptainer build fungiforge.sif env/fungiforge.def       # or: singularity build …

# 3. stage the reference databases once (hours; see §6)
export FUNGIFORGE_DB=/scratch/$USER/fungiforge_db
bash bin/fetch_references.sh

# 4. site config + Lmod module (the share/ wrappers)
cp share/fungiforge-env.sh.example ~/.fungiforge-env.sh   # edit paths / partition / account
module use /path/to/fungiforge/share/modulefiles
module load fungiforge

# 5. run (uses hpc_slurm,apptainer by default)
fungiforge-run --samplesheet samples.csv --run_interproscan true
```

The `hpc_slurm` profile (`conf/hpc_slurm.config`) sets the SLURM executor with time-based
partition routing (`short` for ≤4 h tasks, `long` otherwise) unless you override
`--slurm_partition` / `--slurm_account`; it caps concurrency at 100 queued jobs and rate-limits
submission. The **`share/` Lmod wrappers** — `fungiforge-run`, `fungiforge-env.sh.example`,
and `modulefiles/fungiforge/0.1.0.lua` — put `FUNGIFORGE_HOME`, `FUNGIFORGE_DB`, the work dir,
the profile, and the Apptainer bind paths (`APPTAINER_BINDPATH` must include the DB and work
dirs) in one sourced site file so a run collapses to `fungiforge-run --samplesheet samples.csv`.
Launch from a **writable** directory (Nextflow writes `.nextflow/` and `work/`) and put `work/`
on scratch.

## 4.3 How FungiForge uses containers and conda envs

FungiForge deliberately uses **many small, isolated per-stage environments** rather than one
monolith, because the stage toolchains conflict (a Flye/Medaka assembly stack, a Perl-heavy
Funannotate stack, an antiSMASH stack, RepeatModeler's dfam/tetools stack). Each stage's
`conf/base.config` label points at **both** a conda `env/*.yml` and a container image; the
active packaging profile selects which. Two container details are load-bearing and documented
in the config:

- The **`--data_dir` bind is explicit.** Nextflow auto-mounts the project and work dirs but
  **not** `--data_dir`, so under Docker the profile injects `-v "$data_dir":"$data_dir"` (and
  under Apptainer/Singularity `-B "$data_dir"`) — otherwise the reference DBs are invisible to
  the containers. The path is quoted to survive spaces (e.g. a path containing spaces).
- The **antiSMASH image's entrypoint is cleared** (`--entrypoint=""`) because
  `antismash/standalone` sets `ENTRYPOINT ["antismash"]`, which would hijack Nextflow's bash
  launch; and the Docker profile passes `-e HOME=/tmp -e USER=fungiforge` so tools that call
  `getpass.getuser()`/`expanduser()` (Medaka/torch, Funannotate) do not crash on a host UID
  absent from the container's `/etc/passwd`.

You never activate any of these by hand — Nextflow provisions the right env/container per
process automatically.

**Every image is pinned.** `conf/base.config` names a versioned tag for each public image
(`staphb/flye:2.9.6`, `staphb/spades:4.3.0`, `staphb/medaka:2.2.2`, `staphb/polypolish:0.6.1-bwa`,
`ezlabgva/busco:v5.7.1_cv1`, `dfam/tetools:2.00`) and records, beside each, the manifest digest the
tag resolved to when it was validated; funannotate has no versioned tag for the validated build, so it
is pinned by digest (`nextgenusfs/funannotate@sha256:2eadfeb4…`). The fungiforge base image is built
from `env/base.linux-64.lock`, an explicit URL+md5 lock of the whole solved environment (exported
from the validated image with `bin/lock_env.sh`); `env/base.yml` carries the same versions as exact
`==` pins for reading. Rebuilding from the lock reproduces the conda layer package-for-package
(verified: two independent builds gave identical 453-package lists). `provenance.json` (§11.3)
records what actually ran.

# 5. The `fungiforge` CLI and how to run

## 5.1 The console command

`pip install -e .` registers a small stdlib-only console command, **`fungiforge`** (defined in
`fungiforge/cli.py`, entry point `fungiforge.cli:main`), with four subcommands:

- `fungiforge run …` — thin wrapper over `nextflow run main.nf` (adds `--samplesheet`,
  `--data_dir`, `--outdir`, `-profile`, `-resume`; extra args pass through).
- `fungiforge samplesheet …` — build a sample sheet by pairing ONT and/or Illumina
  read globs by sample id (§7).
- `fungiforge fetch-refs --data_dir DIR` — print the reference-DB fetch command (§6).
- `fungiforge version`.

The raw forms always work unchanged (`nextflow run main.nf …`).

## 5.2 The `nextflow run` invocation and composable profiles

A run composes **one execution profile** and **one packaging profile**:

```bash
# single machine, Docker (on arm64 hosts: linux/amd64 emulated)
nextflow run main.nf -profile local,docker \
    --samplesheet samples.csv --data_dir /path/to/fungiforge_db

# SLURM cluster, native (no emulation)
nextflow run main.nf -profile hpc_slurm,singularity \
    --samplesheet samples.csv --data_dir /path/to/fungiforge_db

# equivalently, via the CLI front door
fungiforge run -profile local,docker \
    --samplesheet samples.csv --data_dir /path/to/fungiforge_db
```

- **Execution:** `local` (this single machine) or `hpc_slurm` (SLURM cluster).
- **Packaging:** `conda` \| `mamba` \| `docker` \| `singularity` \| `apptainer`.
- **`test`** is a third profile that points at the bundled synthetic fixture and is composed
  as `-profile test,docker` (usually with `-stub-run`) to validate the DAG in minutes.

`main.nf` enforces **fail-loud guards** at launch: it errors if `--samplesheet` is missing,
warns if `--data_dir` is unset (identification / resistance / BGC / novelty need the reference
DBs), and errors if `--basecall` is set without `--dorado_model`.

## 5.3 Quickstart (test fixture)

```bash
# validate all 15 stages wire end-to-end in minutes (machinery check, not data quality)
nextflow run main.nf -profile test,docker -stub-run
```

The `test` profile (`conf/test.config`) supplies a subsampled sample sheet, a mini reference
DB, `busco_lineage = saccharomycetes_odb10`, and writes to `results_test/`. Use it to confirm
an install or a config change.

## 5.4 The test suite

| Layer | What | How |
|---|---|---|
| Unit | every helper in `bin/*.py` and the CLI, on synthetic inputs (Kraken2 report fixture, GenBank records built in the test, synthetic proteins / promoters); the scripts are driven as subprocesses the way the pipeline drives them, with subprocess coverage of `bin/` + `fungiforge/` (≥ 85 % enforced) | `pip install -e ".[dev]"`, then `COVERAGE_PROCESS_START=.coveragerc coverage run -m pytest test/unit && coverage combine && coverage report` |
| Static | `nextflow lint` at zero warnings (strict syntax); shellcheck at warning level on every shell script | `nextflow lint .`; `shellcheck -S warning -x bin/*.sh share/bin/*` |
| DAG | the 15 stages on the three-mode fixture under `-stub-run` (routing, joins, DB_CHECK gate, provenance) | `nextflow run main.nf -profile local,test -stub-run` |
| nf-test | `parse_row`/`resolve_path` routing cases (ONT-only, hybrid, Illumina-only, unpaired refused), one stub test per module (every declared output emitted once, named after the sample), the whole DAG with and without the skip flags | `nf-test test` (`nf-test.config`, cases under `test/nf-test/`) |
| Containers | the fixture through the real images with Docker: DB_CHECK and READ_QC complete in-container, versions are recorded, and the assemblers' failure on the synthetic reads is recorded by the status contract, not hidden | `bash test/run_docker_fixture.sh` |
| Real data | CEA10 hybrid and DF-005 Illumina-only compared with `test/expected/*.json` (`fungiforge validate`), assemblies benchmarked against A1163 / NRRL 3357 (`bin/benchmark_assembly.py`), and cyp51A TR34/L98H controls whose genotype is first confirmed from their reads (`bin/control_genotype.py`, `bin/fetch_controls.sh`) | `docs/validation.md` |

The first four layers run on every push (`.github/workflows/ci.yml`); the last two need the
image set and databases and run locally or on the development deployment.

# 6. Reference databases

FungiForge follows the forge convention: **fetch and checksum every database up front, never
mid-run.** `bin/fetch_references.sh` stages them all into one `--data_dir` root using
**aria2c** (16 parallel connections, `--continue` to resume a broken transfer, infinite
retries; `curl -C -` fallback). Each database is guarded by a `.done` sentinel and recorded in
`MANIFEST.tsv`; **one database failing does not abort the others**, and version-sensitive URLs
(UNITE, Kraken2, sourmash-fungi, RVDB) resolve the current release at run time and log-and-skip
rather than guess on a 404. Databases total roughly **150–250 GB**, so they are staged on an
external/scratch drive.

## 6.1 The databases

| Database | Provides | Stage(s) | `fetch_references.sh` step |
|:----------|:----------------------------|:----------|:--------------------------|
| **Container images** | all stage images (Flye, Medaka, antiSMASH 8, Funannotate, dfam/tetools, BUSCO) | all | `images` |
| **antiSMASH** databases (~9 GB) | fungiSMASH BGC detection | 11 | `antismash` (via the antiSMASH 8 image) |
| **Funannotate** DB (~30–50 GB: Pfam, dbCAN, MEROPS, InterPro, BUSCO) | eukaryotic annotation | 07 | `funannotate` (`funannotate setup -i all`) |
| **eggNOG** 5.0.2 (~12 GB compressed: `eggnog.db`, `eggnog_proteins.dmnd`, taxa) | eggNOG-mapper orthology annotation | 07b | `eggnog` (direct download) |
| **InterProScan** data release 5.78-109.0 (6.9 GB; matches the image tag) | InterProScan domains / GO | 07c | `interproscan` (on request) |
| **dbCAN** family HMMs (V14), **PHI-base** FASTA, **EffectorP 3** (repository with WEKA) | CAZymes, virulence, effectors | 13 | `dbcan`, `phibase`, `effectorp` |
| **geNomad** database v1.9 (0.8 GB) | viruses / plasmids / proviruses | 10 | `genomad` |
| **Reference genome set** (NCBI Datasets, one reference genome per species of the genera in `novelty_genera.txt`; tens of GB) | genome-ANI novelty | 12 | `genomes` (on request) |
| **BUSCO / compleasm** `fungi_odb10` + the order/class lineages of `busco_lineages.tsv` | assembly completeness (gate, then species-aware) | 05, 08b | `busco` (`BUSCO_LINEAGES` overrides the set) |
| **UNITE** general FASTA (Fungi v10.0, 2025) | ITS species identification | 08 | `unite` |
| **Kraken2** (PlusPF-8 GB) | decontamination | 04 | `kraken2` |
| **sourmash/RefSeq-fungi** (GenBank fungi k=31) | genome-ANI ID & novelty | 08, 12 | `refseq_fungi` |
| **FungAMR** tables + built reference proteins + derived panel | antifungal-resistance calling (panel rows with evidence tiers are derived from the table at run time when `fungamr_panel.tsv` is absent) | 09 | `fungamr` (+ `build_fungamr_refs.py`) |
| **RVDB-prot** (v31.0, 786 k proteins) | mycovirus / endogenous-viral-element screen | 10 | `rvdb` |
| **Type-material marker sets** (NCBI "sequence from type": CaM, BenA, TEF1, RPB2, LSU; optional) | secondary-locus identification | 08 | `markers` (`fetch_marker_refs.py`; `NCBI_API_KEY` speeds it up) |
| **PubMLST fungal schemes** (optional) | MLST | 08 | `mlst` (`fetch_mlst_schemes.py`) |

**Checked before every run.** The first task of a run, `DB_CHECK`, verifies that every database the
enabled stages need is present and complete (the `.done` marker `fetch_references.sh` writes plus the
key files each stage opens: UNITE FASTA, Kraken2 `*.k2d`, BUSCO lineage, funannotate `Pfam-A.hmm`,
antiSMASH `clusterblast/`, FungAMR tables + reference proteins; `refseq_fungi` only with `--genome_id`;
Kraken2/antiSMASH only when their stages are on). Anything missing is listed once and the run stops
in seconds instead of hours later; `--allow_missing_db` turns that into a warning, and the
dependent stages then record `skipped`. `--funannotate_db` / `--antismash_db` overrides are honoured.
The findings are written to `pipeline_info/db_manifest.json` (with each database's `MANIFEST.tsv`
source line) and embedded in `provenance.json`.

## 6.2 Staging them

```bash
export FUNGIFORGE_DB="/path/to/fungiforge_db"
bash bin/fetch_references.sh                       # all steps (many hours; resumable)
bash bin/fetch_references.sh images antismash fungamr   # only selected steps
```

Re-running is safe and resumes where it stopped. The FungAMR step is special: FungAMR ships a
**mutation table with no sequences**, so after downloading the tables the script calls
`build_fungamr_refs.py`, which reads every `ref_seq_uniprot_accession`, fetches those proteins
from **UniProt** (batched, with recursive batch-splitting so one obsolete accession never sinks
its neighbours), and writes `fungamr/reference_proteins.faa` with `>GENE__ACC__Species` headers
— the sequences Stage 09 aligns against. **MARDy is not fetched** (its REST API is API-key-gated
and its content is largely subsumed by FungAMR).

> **Note on refreshing version-sensitive URLs.** The script embeds the current default URL and
> the recipe to resolve a newer one for UNITE (PlutoF DOI), Kraken2 (genome-idx S3),
> sourmash-fungi (UC Davis prepared DB) and RVDB (Institut Pasteur). If a URL 404s, that step
> logs `FAILED` in the manifest and continues; resolve the new link and re-run just that step.

# 7. Inputs — the sample sheet

The only per-run input is a **CSV sample sheet** (`--samplesheet`), with this exact schema
(parsed in `main.nf`):

```text
sample,ont_fastq,illumina_r1,illumina_r2,rna_r1,rna_r2,compartment,facility,season
```

The two RNA columns are optional and may be omitted entirely; a sheet written before they existed
is still valid.

- **`sample`** — unique isolate id (becomes the output directory and every filename).
- **`ont_fastq`** — the ONT FASTQ. **Optional** — required only when the isolate has no Illumina
  reads. When `--basecall` is set this column is instead the **pod5 directory** Dorado basecalls.
- **`illumina_r1`, `illumina_r2`** — paired short reads. **Optional** individually, but they come
  as a pair (one without the other is an error). With ONT present → the Stage 03b hybrid polish
  runs (under `--hybrid auto`). Without ONT → the isolate is **Illumina-only** and is assembled by
  SPAdes (Stage 02b). Missing reads are substituted with the `assets/NO_ONT` / `NO_R1` / `NO_R2`
  sentinels so the channel plumbing stays well-typed.
- **`rna_r1`, `rna_r2`** — paired RNA-seq for this isolate, used as gene-prediction evidence
  (Stage 07a). **Optional**, but paired like the Illumina columns. Without them prediction stays
  *ab initio*; with them the models are trained on transcripts. Missing values are substituted with
  the `assets/NO_RNA1` / `NO_RNA2` sentinels, and `meta.has_rna` records which applies.
- **`compartment`, `facility`, `season`** — the One Health metadata carried untouched to the
  master table and consumed by Layer 2 (e.g. `AIR` / `HUMAN` / `SURFACE`; a facility label; a
  season). Default to `NA` when absent.

**Each row must provide either `ont_fastq` or the `illumina_r1`+`illumina_r2` pair (or all three).**
The mode is auto-detected per row (`meta.assembly_mode`): ONT present → `longread`
(ONT-only or hybrid); ONT absent → `shortread` (Illumina-only).

A generic example mixing all three modes:

```text
sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season
AF_AIR_001,reads/AF_AIR_001.ont.fastq.gz,reads/AF_AIR_001_R1.fastq.gz,reads/AF_AIR_001_R2.fastq.gz,AIR,Abattoir_1,Dry
CA_HUM_014,reads/CA_HUM_014.ont.fastq.gz,,,HUMAN,Clinic_2,Wet
CA_SUR_009,,reads/CA_SUR_009_R1.fastq.gz,reads/CA_SUR_009_R2.fastq.gz,SURFACE,Clinic_2,Wet
```

Here `AF_AIR_001` is **hybrid** (ONT + Illumina → hybrid-polished, high-confidence resistance),
`CA_HUM_014` is **ONT-only** (provisional resistance), and `CA_SUR_009` is **Illumina-only**
(SPAdes assembly, high-confidence resistance).

**Building one automatically.** `fungiforge samplesheet` pairs ONT and/or Illumina globs
by sample id and stamps the metadata (omit `--ont` entirely for an Illumina-only cohort):

```bash
fungiforge samplesheet \
    --ont 'reads/*.ont.fastq.gz' \
    --illumina-r1 'reads/*_R1*.fastq.gz' --illumina-r2 'reads/*_R2*.fastq.gz' \
    --compartment AIR --facility Abattoir_1 --season Dry \
    -o samples.csv
```

It derives the id by stripping the read suffix, matches R1/R2 to the ONT id (dropping `_R1`/`_1`
markers), and leaves the Illumina columns blank for any sample without short reads. Edit the CSV
afterwards to correct per-sample compartment/facility/season if they differ.

# 8. Configuration & parameters

Set parameters with `--key value` on the command line (there is no params-file convention baked
in, but `-params-file x.yaml` works as usual for Nextflow). Compute caps and stage toggles are
the knobs you will touch most.

**Key parameters (real defaults from `nextflow.config`).**

| Param | Default | Meaning |
|---|---|---|
| `--samplesheet` | `null` | **required** — the CSV sample sheet |
| `--data_dir` | `null` | reference-database root (§6); warns if unset |
| `--outdir` | `results` | results directory |
| `--hybrid` | `auto` | ONT isolates only: `auto` (polish with Illumina if present) \| `on` \| `off` |
| `--assembler` | `flye` | long-read assembler (ONT/hybrid): `flye` \| `canu` \| `raven` |
| `--sr_assembler` | `spades` | short-read assembler (Illumina-only): `spades` \| `megahit` |
| `--purge_dups` | `true` | collapse heterozygous haplotigs after long-read assembly (skipped for short-read) |
| `--min_contig_len` | `500` | drop contigs shorter than this, or with fewer than 4 distinct bases, at short-read assembly and before annotation (`funannotate sort --minlen`); short-read assemblers leave many tiny fragments and occasional homopolymer stubs that funannotate rejects |
| `--ont_min_qual` | `10` | chopper Q filter (Stage 01) |
| `--ont_min_len` | `1000` | chopper length filter (Stage 01) |
| `--basecall` | `false` | run Dorado on a pod5 dir (ont_fastq column = pod5 dir) |
| `--dorado_model` | `r1041_e82_400bps_sup_v5.2.0` | Dorado model; also Medaka's model |
| `--dorado_duplex` | `false` | Dorado duplex basecalling |
| `--busco_lineage` | `auto` | `auto` (order-specific after ID) \| `fungi_odb10` \| `<lineage>` |
| `--genome_id` | `false` | Stage 08 genome-level sourmash gather (slow emulated; ITS is primary) |
| `--skip_rna_train` | `false` | ignore the samplesheet RNA columns and predict *ab initio* |
| `--rna_stranded` | `no` | `funannotate train --stranded`: `no`, `RF`, `FR`, `F`, `R` |
| `--rna_max_intronlen` | `3000` | `funannotate train --max_intronlen` |
| `--funannotate_seed` | `anidulans` | fallback Augustus species when the species map has no entry and when the training helper cannot run |
| `--run_interproscan` | `false` | Stage 07c InterProScan in its own image (heavy; needs `--interproscan_data`) |
| `--interproscan_data` | `null` | InterProScan data directory. Defaults to `<data_dir>/interproscan/data`, the stable symlink `fetch_references.sh interproscan` writes, so only a non-standard location needs this |
| `--skip_eggnog` | `false` | remove Stage 07b eggNOG-mapper |
| `--genemark_dir` | `null` | unpacked GeneMark-ES directory (licensed; bound into the container) |
| `--ploidy` | `auto` | ploidy handling |
| `--ploidy_min_sites` | `10000` | biallelic sites nQuire needs before its ploidy call is trusted (Stage 13) |
| `--genemark_key` | `null` | path to the free-academic GeneMark licence key (with `--genemark_dir`). Its directory is bound into the container, so the key may sit beside the GeneMark directory rather than inside it |
| `--af_panel` | bundled `af_resistance_panel.tsv` | curated resistance panel (Stage 09) |
| `--skip_decontam` | `false` | skip Stage 04 decontamination |
| `--skip_mge` | `false` | skip Stage 10 mobile elements |
| `--skip_bgc` | `false` | skip Stage 11 BGCs |
| `--skip_novelty` | `false` | skip Stage 12 novelty |
| `--skip_extras` | `false` | skip Stage 13 extras |
| `--max_cpus` | `16` | max cores per task on the local machine |
| `--max_memory` | `120.GB` | memory ceiling |
| `--max_time` | `96.h` | time ceiling (fungal annotation is slow, esp. emulated) |
| `--slurm_partition` / `--slurm_account` | `null` | SLURM routing (else time-based `short`/`long`) |

**Stage toggles genuinely skip.** A `--skip_*` flag makes the subworkflow emit an *empty
channel* for that stage rather than running the process emptily — so the stage's container is
never even pulled. Skippable stages are 10 (mobile), 11 (BGC), 12 (novelty), 13 (extras), plus
04 decontam via `--skip_decontam`. The assembly-through-report spine (00–09, 14) is not
skippable.

**Per-process resources.** Resources are declared per stage via labels in `conf/base.config`
(e.g. `assembly` gets up to 12 CPUs / 64 GB / 24 h; `annotate` up to 12 CPUs / 64 GB / 72 h;
`repeat` up to 12 CPUs / 48 GB / 48 h). CPU-bound tools receive `task.cpus`. The default process
policy retries on out-of-memory/time exit codes (137/139/140/143/247) up to twice, scaling
memory and time with the attempt number.

# 9. The stages in detail

Each stage publishes to `results/<sample>/<NN>_<stage>/` and emits a `result.json`. Generic
example commands below show the pipeline invocation that exercises the stage; you normally run
the whole pipeline, not stages individually.

## 9.1 Stage 00 — Basecall (optional)

Runs **Dorado** on a pod5 directory to produce the ONT FASTQ, only when `--basecall` is set
(and `--dorado_model` must be given, or `main.nf` errors). Dorado runs as a **native host
binary**, which keeps basecalling fast on Apple Silicon. Supports duplex via `--dorado_duplex`.

```bash
nextflow run main.nf -profile local,docker --samplesheet pod5_samples.csv \
    --data_dir /path/to/fungiforge_db --basecall --dorado_model r1041_e82_400bps_sup_v5.2.0
```

## 9.2 Stage 01 — Read QC & filter (mode-aware)

For **long-read isolates** (ONT present): **NanoPlot** on raw and filtered reads; **chopper**
filters by quality (`--ont_min_qual`, default Q10) and length (`--ont_min_len`, default 1000 bp);
**Illumina (if present):** fastp. For **Illumina-only isolates**: **fastp** only, and an empty ONT
placeholder is emitted to keep the tuple contract. `readqc.json` records the `platform`
(`ont` \| `hybrid` \| `illumina`) and the filtered ONT FASTQ (long-read) is carried forward.

### Stage 00b — read-level triage (W2.2)

Kraken2 on a subsample of the QC'd reads (`--triage_reads`, default 200,000; Illumina when
present, ONT otherwise) gives an early verdict with the same rules as the contig-level verdict
after assembly. An isolate that is `non_fungal` or `human` at read level is stopped **before
assembly** (gate reason `read_triage:<verdict>`), so a bacterium on a fungal plate costs a few
minutes rather than hours; its master row carries `read_verdict`, the top taxon and the gate
reason. `likely_fungal` and `mixed` proceed. `--read_triage false` disables the stage;
`--force_all` disables the gate. Without a Kraken2 database the verdict is `not_run`.

### Stage 01b — k-mer profile (W2.2)

KMC counts k-mers (`--kmer_k`, default 21) in the trimmed Illumina reads (or the filtered ONT
reads) and GenomeScope2 models the spectrum at ploidy 1 and 2. The stage JSON and master row
report the estimated genome size, heterozygosity, a ploidy hint (diploid only when the ploidy-2
model fits at least as well and heterozygosity is at least 0.5 %), k-mer coverage and read
coverage. Models that do not converge (low coverage, very noisy reads) give null values with a
note and mark the stage `partial`; nothing is invented. Compare `genome_size_est` with the
assembly length in Stage 05: a large gap points at contamination, a collapsed diploid, or a
partial assembly. `--kmer_profile false` disables the stage.

## 9.3 Stage 02 — Assembly (long-read) / Stage 02b — Assembly (short-read)

**Long-read isolates (`assembly_mode == longread`)** are assembled with the chosen `--assembler`:
**Flye `--nano-hq`** (default), Raven, or Canu (`genomeSize=35m`). When `--purge_dups true`
(default), **purge_dups** collapses heterozygous haplotigs via a minimap2 self-alignment; if
purge_dups is unavailable the raw assembly is kept. Organelle contigs are then put back:
purge_dups judges contigs by read depth against the nuclear peak, a mitochondrial genome sits far
above it, and it is duly discarded as collapsed duplication. `protect_organelle.py` restores only
contigs carrying at least two core mitochondrial genes within the length and GC bounds, the same
evidence Stage 04 uses, so genuine haplotigs stay purged. `purge.json` records
`n_organelle_restored`.

**Illumina-only isolates (`assembly_mode == shortread`)** are assembled by **Stage 02b** with
`--sr_assembler`: **SPAdes `--isolate`** (default; built-in read error-correction) or **MEGAHIT**.
purge_dups (a long-read haplotig step) does not apply. Both stages emit the same
`*.assembly.fasta` so the two paths merge transparently downstream.

## 9.4 Stage 03a — ONT polish (Medaka), long-read isolates

**Medaka** (`staphb/medaka`) computes a neural consensus of the draft using the ONT reads and
the `--dorado_model`. If Medaka produces no consensus the draft passes through unchanged. Emits
`*.medaka.fasta`. Runs only for long-read isolates (Illumina-only assemblies skip it — there is no
ONT and no homopolymer error to correct). Medaka and Polypolish live in different containers,
which is why the polish is split into 03a and 03b.

## 9.5 Stage 03b — Short-read polish (Polypolish) / Illumina-only passthrough

For **hybrid isolates** (ONT + Illumina, `--hybrid` `auto`/`on`), this stage polishes the Medaka
consensus with **Polypolish** (`staphb/polypolish`). Polypolish requires **all** alignments per
read (`bwa mem -a`) plus an insert-size filter — the module runs `bwa index`, two `bwa mem -a`
passes, `polypolish filter`, then `polypolish polish`. For **Illumina-only isolates** the SPAdes
assembly is already short-read-derived, so it passes through unchanged. The reported **mode is
honest**:

- `hybrid` — the short-read pass succeeded (`resistance_confidence: high`);
- `illumina_only` — Illumina-only assembly, no ONT (`resistance_confidence: high`; no homopolymer-indel risk);
- `ont_only` — long-read isolate with no Illumina reads (`provisional_ont_only`);
- `hybrid_failed_ont_fallback` — the pass was attempted but failed, so the Medaka consensus is
  kept and the mode says so (`provisional_ont_only`) — never silently claimed as a clean hybrid.

The `polish.json` records the mode, the polisher version, the bp changed, and the
`resistance_confidence` that rides downstream to flag Stage 09 calls.

## 9.6 Stage 04 — Decontamination + organelle split

**Kraken2** (against `--data_dir/kraken2`) classifies every contig; bacterial, archaeal, viral and
human contigs are dropped, the domain composition and a fungal / non-fungal verdict are recorded
(`sample_verdict`, `contam_removed_pct`, `top_taxon`), and a non-fungal isolate is kept whole and
flagged for the gate rather than "decontaminated". The **mitochondrial genome** is then separated
from the assembly (`mito_extract.py`): a contig is mitochondrial when tblastn of the bundled core
mitochondrial proteins (`fungiforge/resources/mito/`, 15 genes from *A. fumigatus*, *S. cerevisiae*
and *C. neoformans*) finds at least two core genes on it, it is at most `--mito_max_len` (250 kb)
and AT-rich (GC < 40 %). Emits `*.nuclear.fasta`, `*.mito.fasta` and `decontam.json` (with the
mitochondrial contigs and rejected candidates).

## 9.6b Stage 04b — Mitochondrial genome

`mito_annotate.py` annotates the separated mitogenome: the 15 core genes by tblastn (best
reference species, identity, coverage, exons and introns from HSP gaps, extra copies), rnl/rns by
blastn, circularity from a terminal repeat, and, from a read subsample mapped to mitogenome +
nuclear genome, the mito/nuclear depth ratio (copy number) and heteroplasmic sites (minor allele
≥ 10 % at ≥ 20× for Illumina, ≥ 20 % for ONT; the pileup uses reads downsampled to
`--organelle_pileup_depth` × (default 300, `0` = every read); `--organelle_reads_check false` skips the
mapping). Emits `organelle.json` and
`*.mito.gff`. tRNAs are not annotated.

## 9.7 Stage 05 — Assembly QC + completeness

**QUAST** contiguity plus **compleasm** (preferred) or **BUSCO** completeness against the
order-specific ODB10 lineage (`fungi_odb10` fallback under `--busco_lineage auto`).
`bin/assembly_qc.py` computes assembly length, contig count, N50, largest contig and the BUSCO
% complete, and sets a **MIMAG-style `qc_pass`** flag (true when the assembly is non-empty and
BUSCO ≥80% or unavailable). Emits `assemblyqc.json`.

### Stage 05b — the gate (W2.1)

After decontamination and assembly QC, each isolate is checked before the expensive fungal
stages: an isolate whose Kraken2 verdict is `non_fungal` or `human`, or whose assembly failed
QC (`qc_pass` false), is **stopped**. It receives a `gate` stage JSON (status `skipped`, reason
`verdict:<verdict>` or `qc_pass:false`), still reaches Stage 14, and so still has a master row —
with `gate = skipped(<reason>)`, `gate:skipped` in `stages_failed`, and `NA` for everything
downstream — while repeat masking, annotation, identification, resistance, mobile elements,
BGC, novelty and extras are not run for it. A bacterial isolate on a fungal plate therefore
costs minutes, not the hours of annotation it used to. `--force_all` sends every isolate through
regardless (for example to inspect a `mixed` sample by eye; `mixed` itself is not stopped).
Isolates whose decontamination did not run (`verdict = not_run`) pass the gate.

## 9.8 Stage 06 — Repeat modelling + soft-masking

**RepeatModeler2** builds a de-novo TE library (`-LTRStruct`), **RepeatMasker** soft-masks the
nuclear genome (`-xsmall`) — both from the `dfam/tetools` image. Emits the soft-masked genome
(for annotation), the **TE library** (for Stage 10), and `repeat.json`.

## 9.9 Stage 07 — Eukaryotic annotation (Funannotate, four processes)

**07a Predict.** `funannotate predict` (GeneMark-ES / Augustus / SNAP / GlimmerHMM combined by
EVM) from the `nextgenusfs/funannotate` image with `FUNANNOTATE_DB=--data_dir/funannotate`.
GeneMark is licensed and not in the image: give `--genemark_dir` (the unpacked GeneMark-ES
directory; the profiles bind it into the container) and `--genemark_key`; without them
prediction runs without GeneMark and `predict.json` says so. Training is **species-aware**: the
stage-08 species call selects the Augustus pre-trained seed species and the funannotate BUSCO set
from `fungiforge/resources/annotation_training.tsv` (`annotation_training.py` uses a mapped
choice only when it is staged, else `anidulans` / `dikarya`). Emits `predict_results/` and the
predicted proteins.

When the samplesheet gives `rna_r1`/`rna_r2`, prediction runs `funannotate train` first
(Trinity, PASA and HISAT2, all present in the image) into the same output directory, which
`funannotate predict` then picks up automatically, followed by `funannotate update` for UTRs and
transcript-corrected models. The updated models replace `predict_results`, so stages 07b, 07c and
07 see transcript-informed proteins rather than stale ones. `predict.json` records `rna_trained`
and `models_from`. Expect this to be slow: PASA's SQLite backend is single-threaded whatever
`--cpus` says. The stage probes for Trinity and PASA rather than assuming them, and degrades to
*ab initio* prediction with a recorded reason if either is missing, if the FASTQ files are empty,
or if training produces no PASA GFF3.

**07b eggNOG-mapper** (`eggnog-mapper` image, eggNOG 5.0.2 data from `fetch_references.sh
eggnog`, fungal taxonomic scope) and **07c InterProScan** (`interpro/interproscan` image with the
data release from `fetch_references.sh interproscan`, found under `--data_dir` without further
configuration; enable with `--run_interproscan true`, off by default because it takes hours per
genome) annotate the predicted proteins in
parallel; each is skipped with an empty table when its data are not staged, and `--skip_eggnog`
removes 07b. Stage 07b treats a database that is present but short (a staging run still downloading
or decompressing) as not staged, and an emapper that fails is recorded `partial` rather than failing
the isolate, so a run started while `fetch_references.sh eggnog` is in flight still completes; add
the table later with `-resume`, which reruns only 07b and annotation.

**07 Annotate.** `funannotate annotate` (Pfam, dbCAN, MEROPS, UniProt, BUSCO, plus `--eggnog` /
`--iprscan` when the tables exist) on a private copy of the prediction. Emits the **proteins
FASTA** and the **GenBank (`.gbk`)** that Stages 09–13 consume, plus `annotate.json` with the
annotation coverage (`annotate_stats.py`: proteins, % with PFAM / InterPro / GO / eggNOG / a named
product, secreted, CAZyme, protease, BUSCO, EC counts). This is the slowest stage under emulation.

## 9.10 Stage 08 — Identification (multi-locus concordance)

`barrnap` finds the rRNA operon; `extract_rrna_region.py` extracts a padded window around it so
**ITSx** stays under the HMMER 100 kb limit (it aborts on whole chromosomes); ITSx carves
ITS1/5.8S/ITS2; **vsearch** classifies ITS against **UNITE** (`--id 0.90`, top hits). The
**secondary loci** are extracted from the assembly by `extract_markers.py`: calmodulin (CaM),
β-tubulin (BenA), TEF1-α and RPB2 are located by `tblastn` with bundled reference proteins
(`fungiforge/resources/markers/`), HSPs on the best contig/strand are chained across introns and
the locus is written with flanks in coding orientation; the LSU D1/D2 region is the first 900 bp
of barrnap's 28S. Each locus is searched (`blastn`, megablast) against a reference set of NCBI
records from **type material** (`--data_dir/markers/<locus>.fasta`, `fetch_references.sh markers`).
**MLST** runs against the fungal PubMLST schemes (`--data_dir/mlst`, `fetch_references.sh mlst`;
*A. fumigatus*, *C. albicans*, *C. glabrata*, *C. tropicalis*, *C. krusei*) through
`mlst_fungal.sh`. Optionally (`--genome_id true`, slow emulated) **sourmash gather** does
genome-level ANI.

`id_classify.py` resolves the call: **high** needs the ITS species *and* at least one agreeing
secondary line (locus or genome match) with none contradicting; **medium** is ITS alone
(flag `secondary_loci_unavailable` when no reference set is staged), secondary-only, ITS at
94–98.5 %, or ties that include the ITS species; **low** is genus-only or a discordance (the call
becomes *Genus* sp. with `discordant:<locus>=<species>` in `flags`). A locus is a *tie* when
reference species within 0.3 % identity of the best hit differ (e.g. *A. flavus* vs *A. oryzae* on
CaM and BenA): it is reported with its candidates and never confirms a species by itself. The
species-aware BUSCO lineage is chosen from `fungiforge/resources/busco_lineages.tsv` (species row,
then genus, else `fungi_odb10`) and written to `*.busco_lineage.txt`. Emits `*.species.txt`,
`*.markers.fasta` (ITS, rRNA window and every extracted locus), `identify.json` (`loci` with
per-locus level/identity/candidates, `concordance`, `flags`, `mlst`, `busco_lineage`).

## 9.10b Stage 08b — Species-aware BUSCO

Stage 05 scores completeness against `fungi_odb10` before the species is known (the QC gate);
stage 08b re-scores the nuclear assembly with the lineage stage 08 chose (e.g. `eurotiales_odb10`
for *Aspergillus*, `saccharomycetes_odb10` for *Candida*). It records `skipped` with the reason when
the lineage is the one stage 05 used, when it is not staged under `--data_dir/busco`
(`fetch_references.sh busco` downloads the common lineages straight from the BUSCO data server into `--data_dir/busco/lineages/`; `BUSCO_LINEAGES` overrides the set), or
when `--busco_lineage` is fixed. Emits `busco_lineage.json` (`lineage`, `busco_complete`, S/D/F/M);
the master row carries `busco_lineage_specific` and `busco_complete_specific` next to the
stage-05 values.

## 9.11 Stage 09 — Antifungal resistance (the centerpiece)

The bespoke core; see §10 for full depth. `af_resistance.py` takes the Funannotate proteins, the
species call, the nuclear assembly and the GBK; selects species-relevant panel rows from the
curated overlay **merged with rows derived from the FungAMR catalogue** (`fungamr_panel.py`:
every published substitution with its evidence tier, 1 = engineered and measured in the same
species … 8 = seen in a resistant isolate without validation); finds the isolate ortholog of each
target gene; reads residues at the panel positions; classifies each as a known resistance
mutation, an `associated_unvalidated` tier-8 change, a novel hotspot variant, wild-type,
loss-of-function, gain-of-function, or over-expression target; and runs the `cyp51a_TR.py`
promoter-TR detector for *A. fumigatus*.

**Read-level genotyping.** A subsample of the QC'd reads (`--genotype_reads` Illumina pairs,
`--genotype_reads_ont` ONT reads) is mapped with minimap2 to the GenBank records
(`gbk_to_fasta.py`, the annotation's coordinates) and `read_genotype.py` translates every read
at each hotspot codon: allele frequency, zygosity (`homozygous` / `heterozygous` / `mixed` /
`insufficient` below `--genotype_min_reads`), agreement with the assembly, known alleles the
assembly lacks (reported as *reads only*), the cyp51A TR site as an insertion or deletion in the
spanning reads (an extra copy the assembly missed becomes a TR call flagged `discordant`), and
locus/genome depth as a copy-number ratio for target and efflux genes. `--read_genotype false`
keeps the assembly-only path.

**Confidence** combines the polish mode, the read support and the evidence tier: `high`
(read-confirmed, or hybrid/Illumina assembly), `provisional_ont_only`, `discordant` (reads
contradict the assembly), `low_evidence` (tier 8), `provisional_minor_allele` (reads-only below
80 %). Emits `resistance.json` (`calls[].read_support`, `calls[].copy_number`, `panel`,
`summary.read_support`, `summary.copy_number_flags`) and the mapped subsample as
`<sample>.reads.bam`.

## 9.12 Stage 10 — Mobile & repeat elements

Four blocks, each optional: (1) the **TE landscape** from RepeatMasker's summary table
(`repeatmasker_tbl.py`: `te_percent`, interspersed-repeat fraction, LTR/Gypsy/Copia, LINE, DNA
transposon, Helitron, unclassified and simple-repeat percentages) plus the RepeatModeler family
counts by superfamily (`te_summary.py`); (2) **geNomad** on nuclear + mitochondrial contigs
(viruses, proviruses, plasmids; `--data_dir/genomad_db`, `fetch_references.sh genomad`;
`--genomad false` disables); (3) the **mycovirus / endogenous-viral-element screen**: DIAMOND
blastx of every contig against RVDB-prot (`--data_dir/rvdb`, database built in the task;
`--rvdb_evalue`, `--rvdb_sensitivity`), loci classified by virus family — mycovirus families,
retroelement-like (LTR retrotransposons, reported apart), other viral — with the DNA-only caveat:
RNA mycoviruses are visible only as reverse-transcribed endogenous copies; (4) **mitochondrial
homing-endonuclease ORFs** (`mito_heg.py`: six-frame ORFs under genetic code 4 searched with the
Pfam LAGLIDADG / GIY-YIG / HNH profiles from the funannotate database). Emits `mobile.json` and
the raw RVDB hits.

## 9.13 Stage 11 — Biosynthetic gene clusters

**antiSMASH 8** in fungal mode (`--taxon fungi`, gene finding off: the Funannotate GenBank carries
the genes) with **KnownClusterBlast** against MIBiG (`--antismash_extra`, default
`--cb-knownclusters`). `bgc_summary.py` counts regions by product type and records each region's
best MIBiG cluster (accession, description, core-gene hits), flagging known **mycotoxin**,
antifungal/antibiotic and other bioactive clusters from `fungiforge/resources/mycotoxin_compounds.tsv`
(`mycotoxin_clusters`, `bioactive_clusters`, `n_bgc_mibig_hits` on the master row). Best-effort:
an antiSMASH failure is recorded, not fatal. Emits `bgc.json` and the region GenBanks, which the
cohort stage clusters into **gene-cluster families** across isolates (`cohort/bgc_families.tsv`:
regions sharing ≥ `--gcf_min_similarity` of their proteins by DIAMOND homology).

## 9.14 Stage 12 — Novelty

**skani** ANI of the nuclear assembly against the staged reference genome set
(`--data_dir/refseq_fungi_genomes`: the reference genome of every species in the genera of
`fungiforge/resources/novelty_genera.txt`, `fetch_references.sh genomes`; aligned-fraction floor
`--novelty_min_af`). `novelty_call.py` then rules: ANI ≥ 95 % → `known_species` (flagged when the
Stage 08 species differs from the nearest genome's); 90–95 % → `candidate_novel_species` unless
ITS and a secondary locus agree at high confidence with that species; < 90 % or no aligned
reference → `candidate_novel_or_unrepresented`; without a genome set the ITS identity alone
decides (< 98.5 % → candidate novel). Emits `novelty.json` with the five nearest genomes.
Skippable with `--skip_novelty`.

## 9.15 Stage 13 — Eukaryote extras

`extras.py` runs five optional blocks on the annotated proteins, the GenBank and the read BAM
of Stage 09; a missing tool or database leaves that block `null` with the reason (stage
`partial`, never failed):

- **CAZymes** — `hmmsearch` against the dbCAN family HMMs (`--data_dir/dbcan`,
  `fetch_references.sh dbcan`) with run_dbcan's filter (E ≤ 1e-15, coverage ≥ 0.35).
- **Secretome and effectors** — SignalP 6 (licensed: build the site image with
  `bin/hpc_install.sh --signalp <tarball>`, Appendix A) on all proteins, then EffectorP 3
  (`fetch_references.sh effectorp`) on the signal-peptide proteins.
- **Virulence** — diamond against PHI-base (`fetch_references.sh phibase`), best hit per protein
  (identity ≥ 40 %, coverage ≥ 50 %), counts by PHI-base phenotype.
- **Mating type** — Pfam profiles from the funannotate database: an alpha-box protein (PF04769)
  = MAT1-1; an HMG-box protein (PF00505) within six genes of APN2/SLA2 in the GenBank = MAT1-2;
  both = homothallic or a heterozygous diploid.
- **Ploidy** — nQuire on the read BAM (diploid / triploid / tetraploid likelihoods;
  `haploid_like` when too few biallelic sites remain), next to the k-mer hint of Stage 01b.
  Density decides, not goodness of fit: below `--ploidy_min_sites` biallelic sites (10,000 by
  default) the models cannot be separated, so the call is reported as haploid/homozygous at low
  confidence and nQuire's own best fit is kept in the JSON under `model_best` for inspection. A
  30 Mb haploid genome yielding 2,491 sites was otherwise reported as tetraploid.

Emits `extras.json` (`mating_type`, `ploidy`, `n_secreted`, `n_effectors`, `n_cazymes`,
`n_phibase_hits`, per-protein tables, `tools`, `skipped`).

## 9.15b Stage 16 — Cohort phylogenomics and clonality (run-level)

Runs once over every isolate that passed the gates (`--skip_cohort` removes it) and writes
`<outdir>/cohort/`: **species clusters** from skani all-vs-all ANI (`ani_matrix.tsv`,
`species_clusters.tsv`; single linkage at `--cohort_ani_cluster`); a **phylogenomic tree** from
the single-copy BUSCO proteins Stage 05 now exports (genes single-copy in ≥ `--cohort_min_frac`
of the isolates, MAFFT, gap-trimmed supermatrix `cohort_alignment.faa` with partitions, IQ-TREE 3
LG+G with 1000 ultrafast bootstraps or FastTree, `cohort.treefile`; needs ≥ `--cohort_min_isolates`);
and **clonality** within each species cluster: assembly-based core SNP distances to the most
complete member (minimap2 asm5 + paftools.js call, SNVs only) in `snp_distances.<cluster>.tsv`
and clonal groups at `--clonal_snp_threshold` in `clonal_groups.tsv`. `cohort.json` summarises
all of it for Layer 2's transmission objective.

## 9.16 Stage 14 — Report + master table

`make_report.py` merges every per-stage `result.json` into a **self-contained per-isolate HTML
report** (Jinja2 template `fungiforge/report_templates/isolate_report.html.j2`: headline call,
resistance table with read support and evidence tiers, reads and assembly, identification with the
per-locus evidence, annotation and biology, secondary metabolism, mobile elements and organelle,
run status) and writes the isolate's **master row** to `results/04_summary/<sample>.master.tsv`.
Missing stages become `NA` — never a hard failure.

## 9.17 Stage 17 — Cohort summary (run-level)

`cohort_report.py` runs once after every isolate: it merges the per-isolate rows into
**`results/04_summary/master_fungi.tsv`** in schema order, **validates** it against
`fungiforge/resources/master_schema.json` (§11.2) and writes **`cohort_report.html`** (composition
by species / compartment / facility, gates and QC, resistance and mycotoxin tallies, metric ranges,
and the stage-16 clusters, tree, clonal groups and BGC families), **`cohort_summary.json`** (the
same numbers, machine-readable) and **MultiQC custom content** under `04_summary/multiqc/`, running
MultiQC itself when it is in the image. A schema violation fails this stage loudly at the end of
the run without touching the per-isolate results.

# 10. The antifungal-resistance module in depth

This is where FungiForge adds the most value, because there is **no ResFinder for fungi**. The
flow below is `bin/af_resistance.py` (calling `bin/cyp51a_TR.py`), driven by the curated panel and
the FungAMR references.

![Stage 09 antifungal-resistance calling flow: species-relevant panel rows drive ortholog detection against species-matched FungAMR references; observed residues at panel hotspots are classified by mechanism; the cyp51A promoter TR detector runs structurally for *A. fumigatus*; the polish mode sets the confidence flag.](figures/resistance_flow.png){width=100%}

**1 — The curated panel (`af_resistance_panel.tsv`).** A tab-separated interpretation overlay,
one row per gene×species, with columns `gene, organism_regex, drug_class, drugs, mechanism,
hotspot_aa, known_mutations, note, source`. It covers the clinically important axes:

- **Azoles:** *A. fumigatus* **cyp51A** (G54/L98H/Y121F/G138C/M220/T289A/G448S…), **cyp51A
  promoter** (TR34/TR46/TR46³/TR53), **hmg1**, cyp51B; *Candida albicans* / *auris* /
  *parapsilosis* / *tropicalis* / *Nakaseomyces glabratus* / *Cryptococcus neoformans* **ERG11**
  (Y132F, K143R, F126L…); the azole efflux regulators **TAC1/MRR1/UPC2/PDR1** (GOF) and pumps
  **CDR1/CDR2/MDR1/AtrF** (over-expression).
- **Echinocandins:** **FKS1** HS1/HS2 (S645P, F641S…) and **FKS2** across *Candida*/*Nakaseomyces*
  and *A. fumigatus* (S678P).
- **Pyrimidines (5-FC):** **FUR1** (R101C) and loss-of-function of **FCY1/FCY2**.
- **Polyenes:** loss-of-function of **ERG3/ERG6/ERG2**.

**2 — Species-relevance filter.** Only rows whose `organism_regex` matches the isolate's species
(or genus, or `spp.`) are searched — so *A. fumigatus* isolates are not mis-tested against
*Candida* ERG11 numbering.

**3 — Ortholog detection and residue mapping.** For each relevant gene the caller pulls the
matching **FungAMR reference from the isolate's own species** (headers are `GENE__ACC__Species`;
species-matched, else genus, else longest), because the panel's hotspot numbering is
species-specific and a wrong-species reference would misnumber every residue. It finds the
isolate ortholog among the Funannotate proteins with a **k-mer prefilter → BLOSUM62 global
pairwise alignment** (≥55% identity), then reads the residue in the isolate frame aligned to each
1-based reference hotspot.

**4 — Mechanism-aware calling.** Each panel row is classified by its `mechanism`:

- **substitution** — report `known_resistance_mutation` (residue matches a known mutant),
  `novel_hotspot_variant` (changed but not a known mutant), or `wild_type`.
- **loss_of_function** — flag a candidate LoF when the ortholog is truncated (<85% of reference
  length).
- **GOF** (regulators) — report gene presence plus any known gain-of-function residues (confirm
  with expression/phenotype).
- **overexpression** (efflux) — report presence only; over-expression needs RNA-seq and is not
  callable from a genome.

**5 — The cyp51A promoter TR detector (`cyp51a_TR.py`).** For *A. fumigatus* with an assembly, the
structural detector locates cyp51A in the GBK, extracts the ~600 bp upstream promoter (strand-aware),
scans for a tandem array (unit lengths 20–60 bp, ≥90% identity between copies, ≥2 copies), and snaps
the unit length to the nearest canonical TR (TR34/TR46/TR53; TR46³ when ≥3 copies of the ~46 bp
unit). A detected TR is emitted as a high/medium-confidence azole-resistance call annotated with its
canonical coding-mutation pairing (TR34+L98H, TR46+Y121F/T289A).

**6 — The confidence flag (provisional vs high).** Every call carries a confidence derived from the
polish mode: **`high`** for a hybrid-polished assembly, **`provisional_ont_only`** for ONT-only.
This is the pipeline's central honesty guarantee — an ONT-only substitution or frameshift call could
be a homopolymer artefact, so it is never presented as confirmed. The promoter TR call (a structural
signal) is high-confidence when hybrid and medium otherwise.

**7 — Graceful degradation.** Without staged FungAMR references the gene is reported as
`no_reference` ("searched but unresolved") rather than guessed, and the structural TR scan still
runs. The output `resistance.json` summarises the resistant drug classes, the count of known
mutations, and whether a reference was available.

> **Interpretation caveat (in the panel header and every report).** These are **genotype** calls
> that flag *known* resistance-associated changes. **MIC / AFST phenotype confirmation is always
> required**, and over-expression/GOF mechanisms need expression data.

# 11. Outputs

## 11.1 Per-isolate outputs

Under `results/<sample>/`:

| Output | Location | Notes |
|:------------|:--------------|:--------------------------------------------|
| Filtered reads | `01_readqc/` | ONT filtered FASTQ + NanoPlot |
| Polished assembly | `03_polish/` | `*.polished.fasta` + `polish.json` (mode) |
| Nuclear / mito | `04_decontam/` | `*.nuclear.fasta`, `*.mito.fasta`, `organelle.json`, `*.mito.gff` |
| Assembly QC | `05_assembly_qc/` | `assemblyqc.json` (N50, BUSCO, `qc_pass`) |
| Annotation | `07_annotate/` | `*.proteins.faa`, `*.gbk` |
| Identification | `08_identify/` | `*.species.txt`, `identify.json` (loci, concordance, MLST), `busco_lineage.json` |
| Resistance | `09_resistance/` | `resistance.json` (calls + confidence + TR) |
| Mobile / BGC / novelty / extras | `10_…`–`13_…` | per-feature `result.json` |
| **Per-isolate report** | `14_report/` | `*.report.html` (self-contained) |
| **Master row** | `14_report/` and `04_summary/` | `*.master.tsv` |
| **Cohort** | `cohort/` | `ani_matrix.tsv`, `species_clusters.tsv`, `cohort.treefile`, `snp_distances.*.tsv`, `clonal_groups.tsv`, `bgc_families.tsv`, `cohort.json` |

## 11.2 The master table schema (the Layer-2 contract)

The schema lives in `fungiforge/resources/master_schema.json` (name, type, allowed values,
description and group for every column) and is enforced by `bin/validate_master.py`, which Stage 17
runs over the merged table. Columns are **append-only**: a newer pipeline adds columns at the end,
never renames or reorders them, so an older consumer keeps working and a newer table still
validates (extra columns are a warning, not an error). `make_report.py` writes the same column set
to `results/04_summary/<sample>.master.tsv`:

```text
sample, compartment, facility, season,
species, species_confidence, id_method,
qc_pass, busco_complete, busco_lineage, assembly_len, n_contigs,
ploidy, mating_type,
resistant_classes, n_known_af_mutations, cyp51A_TR,
novelty, n_bgc, n_mycovirus, te_percent,
polish_mode,
sample_verdict, contam_removed_pct, top_taxon,
stages_failed,
gate,
read_verdict, genome_size_est, heterozygosity_pct, ploidy_hint, coverage,
id_loci_agree, id_flags, mlst_st, busco_lineage_specific, busco_complete_specific,
resistance_read_support, copy_number_flags,
n_proteins, pct_pfam, pct_go, pct_eggnog, pct_interpro, annotation_training,
n_secreted, n_effectors, n_cazymes, n_phibase_hits,
mito_size_kb, mito_core_genes, mito_circular, mito_copy_ratio, mito_heteroplasmic_sites,
te_ltr_pct, n_genomad_virus, n_genomad_plasmid, n_mito_heg,
mycotoxin_clusters, bioactive_clusters, n_bgc_mibig_hits
```

### Stage status contract

Every per-stage JSON (`results/<sample>/NN_<stage>/<sample>.<stage>.json`) carries a `status`
and a `tools` block, written by `bin/ff_status.sh` / `bin/ff_status.py` at the end of the task:

| `status` | Meaning |
|---|---|
| `ok` | every tool that ran exited 0 |
| `partial` | an optional step failed and the documented fallback was used (e.g. Polypolish failed, the medaka assembly was carried forward, resistance confidence dropped to provisional) |
| `failed` | a required tool failed. The task exits non-zero and the run stops after pending tasks, except for the best-effort stages (BGC, extras), which record `failed` and let the run continue |
| `skipped` | the stage was deliberately not run |

`tools` records each tool's exit code, wall time and (via `ff_version`) the version string the tool
reports; `versions` lists every version recorded in the task, including tools that were probed but not
run; `skipped_tools` records steps not run and why (a missing database, not applicable to this isolate). The master table summarises this per isolate in
`stages_failed` (`stage:status;…` or `none`), so a row with results is never mistaken for a clean run.
A stage that ran no tool at all because its database was absent records `skipped` and appears there
too: an empty eggNOG or InterProScan table is a gap in the row, not a success.
The assembly and medaka steps, which have no scientific JSON of their own, emit small status JSONs
(`<sample>.assemble.json`, `<sample>.medaka.json`) for the same reason. `|| true` is not used
anywhere in the task scripts. The `gate` column (Stage 05b) is `pass`, or `skipped(<reason>)`
for an isolate the gate stopped after assembly QC.

One row per isolate; missing values are `NA` (or `none` for `resistant_classes`). This schema is
deliberately organism-agnostic and mirrors the bacterial study's master table so the two One
Health analyses are directly comparable.

## 11.3 Provenance

Every run ends with a run-level `PROVENANCE` task that writes `pipeline_info/provenance.json`
(schema `fungiforge/resources/provenance.schema.json`, checked at write time):

| Block | Content | Source |
|---|---|---|
| `pipeline` | name, version, git commit and dirty flag, revision/repository when run from GitHub, script hash | `main.nf` (head node) |
| `nextflow`, `run` | Nextflow version/build, session id, run name, start time, full command line, launch/work dirs, profile, container engine, stub/resume flags, user, host | Nextflow `workflow` metadata |
| `params` | every effective parameter | Nextflow |
| `containers` | process label → container reference that was configured | Nextflow |
| `images` | for each reference, what actually ran: a site `.sif` or cached `.img` with size and **sha256**, or a Docker image id + repo digest; `used` says whether any task ran with it (the trace records each task's container, so a base entry overridden by the site config is listed but not counted) | `workflow.onComplete` on the head node (`make_provenance.py images`) |
| `databases` | the `DB_CHECK` manifest: each database's key files, `.done` timestamp and `MANIFEST.tsv` source | `bin/check_databases.py` |
| `samples` | per isolate, per stage: `status`, `tools` (exit, seconds, version), `skipped_tools`, `versions`, `note` | every stage JSON |
| `tool_versions`, `version_conflicts`, `stage_summary` | one consolidated version per tool across the run (a tool reporting two different versions is listed under `version_conflicts` — it should be empty), and ok/partial/failed/skipped counts per stage | aggregated |

Versions are what each tool prints (e.g. `vsearch v2.29.1_linux_x86_64`); the authoritative package
versions of the base image are in `env/base.linux-64.lock`. Hashing the ~30 GB image set takes one to
two minutes at the end of the run; `--provenance_hash_images false` records path and size only. On a
shared install set `params.image_cache_dir` (site config) to the pre-pulled image cache so cached
images resolve; the launcher's `NXF_APPTAINER_CACHEDIR` is used otherwise. Nextflow also writes
`pipeline_info/` (execution report, timeline, trace).

# 12. Layer 2 — the One Health comparative analysis

`analysis/` is a set of **R objective scripts** that consume only the master rows — the same
decoupling the bacterial study uses. Run after the pipeline has produced master TSVs:

```bash
analysis/run_downstream.sh results/ [sample_metadata.csv]
```

Flow: `01_merge_metadata.R` gathers every `04_summary/*.master.tsv` (falling back to
`*/14_report/*.master.tsv`) and optionally left-joins richer sample metadata (air volume,
health-seeking, …) by sample id; `02_recover_metadata.R` standardises compartment/env and derives
resistance/novelty/BGC logicals into `output/analysis_dataset.csv`, which every objective reads. The
eight objectives:

| Script | Objective |
|---|---|
| `03_objective1_composition.R` | community composition & ID (richness, compartment PERMANOVA, priority pathogens) |
| `04_objective2_resistome.R` | antifungal resistome by compartment (azole/echinocandin prevalence, cyp51A TR carriage, Fisher tests) |
| `05_objective3_mobile.R` | TE load & mycovirus carriage by compartment |
| `06_objective4_transmission.R` | cross-compartment clonality (candidate transmission; SNP/phylo confirm) |
| `07_objective5_risk.R` | mixed-effects risk drivers (facility random effect) of resistance/priority carriage |
| `08_objective6_virulence.R` | virulence & resistance×virulence convergence (e.g. azole-R *A. fumigatus* from air) |
| `09_objective7_connectivity.R` | shared BGC / mobile-unit connectivity across compartments (BiG-SCAPE + skani) |
| `10_objective8_novelty.R` | candidate novel species + novel BGCs / production potential |

Conventions match the bacterial study: a shared `_theme.R` with the fixed AIR/HUMAN/SURFACE
palette, `set.seed(20260729)`, self-resolving `output/objectiveN` + `figures/objectiveN`, and
`FigureN_*` PNG+PDF at 300 dpi. The outputs are the substrate for the manuscript and seminar deck
(built with the `scideck` engine).

> **Status.** The objective scripts are scaffolds: several (connectivity, some risk models) emit
> summary tables and note where the full analysis (BiG-SCAPE GCF networks, SNP/phylogenetic
> confirmation) plugs in. They run against real master rows without modification.

# 13. Interpreting results

**Machinery-proven vs needs-real-data.** On the `test` fixture and `-stub-run` the DAG *runs* and
the JSON contracts *fire*, but the numbers are not biology. Treat as machinery-proven only until a
real isolate is processed.

**Read every resistance call through its confidence flag.** A `provisional_ont_only` substitution
or frameshift may be a homopolymer artefact — confirm with a hybrid polish (add Illumina reads and
re-run) before reporting it. The **cyp51A promoter TR** is structural and more robust to ONT error,
but still confirm the paired coding mutation. **Every resistance call is a genotype** — MIC/AFST is
required for a clinical statement.

**Read species calls through confidence and method.** `high`/`GCPSR` (ITS + genome ANI concordant)
is trustworthy; `medium`/`low` (ITS below species threshold, or genus-only) needs multi-locus
confirmation. A `candidate_novel_species` verdict is a hypothesis for polyphasic follow-up, not a
described species.

**Know the current placeholders.** Decontamination keeps all contigs (BlobTools/organelle
extraction pending), efflux over-expression and ploidy are reported as presence/`NA` (need
RNA-seq/read-based tools), and MARDy is not fetched. These are documented honestly in the code and
this manual; do not over-read them.

# 14. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `Missing --samplesheet` at launch | provide `--samplesheet s.csv` with the seven columns |
| `--basecall set but --dorado_model missing` | pass `--dorado_model <model>` |
| ID / resistance / BGC empty, warning about `--data_dir` | set `--data_dir /path/to/fungiforge_db` and stage the DBs (§6) |
| Reference DBs invisible inside a container | the profile binds `--data_dir` automatically; ensure the path is correct and (Apptainer) in `APPTAINER_BINDPATH` |
| Stage 07/11 take many hours locally | linux/amd64 emulation on arm64 hosts — run on a cluster, or accept the wall time; keep `--run_interproscan false` locally |
| antiSMASH "command not found" / container exits immediately | the `bgc` label clears the image entrypoint (`--entrypoint=""`); keep that in `conf/base.config` |
| Medaka/Funannotate crash on `getpass.getuser()` | the Docker profile sets `-e HOME=/tmp -e USER=fungiforge`; keep it |
| ITSx aborts "over comparison pipeline limit" | it was handed a whole chromosome — `extract_rrna_region.py` windows the rRNA operon first; ensure barrnap ran |
| Resistance genes all `no_reference` | FungAMR references not staged — run `bin/fetch_references.sh fungamr` (builds `reference_proteins.faa` from UniProt) |
| A UNITE/Kraken2/RVDB fetch step logs `FAILED` | version-sensitive URL 404'd — resolve the current release (recipe in `fetch_references.sh`) and re-run just that step |
| Resistance calls all `provisional_ont_only` | the isolate is ONT-only — add Illumina reads (populate `illumina_r1/r2`) for a hybrid polish, or run it Illumina-only, for `high` confidence |
| Illumina-only isolate not assembling / wrong assembler | check the row has both `illumina_r1`+`illumina_r2` and no `ont_fastq`; switch `--sr_assembler spades\|megahit` if SPAdes runs out of memory (MEGAHIT is leaner) |
| A stage's container gets pulled despite a skip | ensure you used the `--skip_*` flag (it emits an empty channel, avoiding the pull) |
| conda env won't solve on Apple Silicon | `CONDA_SUBDIR=osx-64 mamba env create -f env/<name>.yml` (some builds are linux/osx-64 only) |

# 15. Reproducibility & provenance

Every image is pinned to a versioned tag or digest in `conf/base.config` (§4.3); the base image
is built from the explicit conda lock `env/base.linux-64.lock`; the manifest
(`manifest.nextflowVersion = '>=23.10'`) pins the engine floor. Databases are verified before every
run (`DB_CHECK`, §6.1) and recorded in `--data_dir/MANIFEST.tsv` with their source URL; FungAMR
reference building records `reference_build.json`. Every run writes `pipeline_info/provenance.json`
(§11.3): commit, Nextflow version, effective parameters, the sha256 or digest of every image that
ran, the database manifest, and every tool version per sample and stage. Cite the image digests and
database releases recorded there. The One Health analysis pins `set.seed(20260729)`.

# 16. Extending FungiForge

- **New organism / study:** change the sample sheet metadata and `--busco_lineage`; no code edits.
- **New resistance knowledge:** add rows to `af_resistance_panel.tsv` (gene×species×mechanism×
  hotspots×known mutations) and rebuild FungAMR references — the caller picks them up.
- **New assembler / polish:** `--assembler {flye,canu,raven}`, `--hybrid {auto,on,off}`.
- **Roadmap (honest):** BlobTools decontam + organelle extraction (Stage 04), read-based ploidy and
  full secretome (Stage 13), BiG-SCAPE GCF networks and SNP/phylo transmission confirmation (Layer 2),
  and a real-SRA end-to-end test.

# 17. Glossary

**ONT** — Oxford Nanopore long-read sequencing. **Hybrid polish** — correcting an ONT assembly with
Illumina short reads (Polypolish). **ITS** — internal transcribed spacer, the fungal barcode.
**UNITE** — curated fungal ITS reference (species hypotheses, `SH…FU`). **GCPSR** — Genealogical
Concordance Phylogenetic Species Recognition. **ANI** — average nucleotide identity. **cyp51A** —
the azole target (14-α sterol demethylase / ERG11). **TR34/TR46** — tandem repeats in the cyp51A
promoter conferring environmental azole resistance. **FKS1/2** — echinocandin target (glucan
synthase), hotspots HS1/HS2. **GOF / LoF** — gain-/loss-of-function. **BGC** — biosynthetic gene
cluster (NRPS/PKS/terpene/RiPP). **fungiSMASH** — antiSMASH fungal mode. **BUSCO / compleasm** —
single-copy-ortholog completeness. **MIMAG** — genome-quality reporting standard. **Master table** —
the fixed-schema per-isolate TSV that is the Layer-1 → Layer-2 handoff.

# 18. Citation

Cite via `CITATION.cff` in the repository root (Bright Adu; *FungiForge: reproducible fungal
genomics from ONT (+ hybrid Illumina)*, v0.1.0, MIT, `github.com/aduton1000/fungiforge`).
FungiForge is the fungal sibling of the forge family — CaptureForge (design), CallForge (targeted
germline analysis), MethylForge (methylomes), and FungiForge (fungal isolate genomics + One Health
comparison).
