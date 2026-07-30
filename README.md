# FungiForge

Reproducible **fungal genomics from Oxford Nanopore** (with an optional Illumina hybrid branch): assembly → eukaryotic annotation → ITS/genome identification → **antifungal-resistance** calling → mobile/TE/mycovirus → **biosynthetic gene clusters** → novelty → One Health comparative analysis.

The fungal sibling of the *forge* family (captureforge / callforge / methylforge). Organism-agnostic by design: **change inputs and config, never code.**

## Two layers

1. **Layer 1 — `fungiforge` (Nextflow DSL2)** — per-isolate genomics, 15 stages (`modules/stage00..14`):
   basecall QC → read QC → assembly → polish (ONT + optional hybrid) → decontam/organelle → assembly QC/BUSCO → repeat-mask → Funannotate → identify → **antifungal resistance** → mobile/TE/mycovirus → **fungiSMASH BGCs** → novelty → eukaryote extras → per-isolate report + `master_fungi.tsv`.
2. **Layer 2 — `analysis/`** — objective scripts (mirroring the bacterial `r-analysis/`) that consume the master table and produce the compartment (air/human/surface) comparative story + figures + manuscript.

## Quickstart (this Apple-Silicon Mac)

```bash
# 1. one-time: fetch reference databases onto the external SSD (hours)
export FUNGIFORGE_DB="/path/to/fungiforge_db"
bash bin/fetch_references.sh                 # aria2c: multi-connection, resumable

# 2. build a samplesheet (ONT required; Illumina optional -> hybrid polishing)
pip install -e .
fungiforge samplesheet --ont 'reads/*.ont.fastq.gz' \
    --illumina-r1 'reads/*_R1*.fastq.gz' --illumina-r2 'reads/*_R2*.fastq.gz' \
    --compartment AIR --facility Abattoir_1 --season Dry -o samples.csv

# 3. run
nextflow run main.nf -profile local,docker \
    --samplesheet samples.csv --data_dir "$FUNGIFORGE_DB"
```

On the HPC (once it's back): `-profile hpc_slurm,singularity` (native, no emulation).

## Profiles

Compose one from each group:
- **execution:** `local` | `hpc_slurm`
- **packaging:** `conda` | `mamba` | `docker` | `singularity` | `apptainer`
- **test:** `test` — subsampled fixture; validate the DAG with `-profile test -stub-run`.

## Notes & best practices baked in

- **ONT resistance-call caveat:** homopolymer indels create false frameshifts exactly where antifungal-resistance point-mutations live. Medaka (+ Illumina Polypolish/POLCA when present) is applied, and resistance calls are flagged **provisional until hybrid-polished**.
- **Antifungal resistance** uses a curated **FungAMR/MARDy** allele+mutation panel plus a dedicated *A. fumigatus* **cyp51A TR34/TR46** promoter-repeat detector (structural, not SNP).
- **Identification** uses ITS/LSU (UNITE) + genome ANI (sourmash/skani) with multi-locus **GCPSR** concordance — there is no GTDB for fungi.
- **Storage:** databases (~150–250 GB) and Nextflow `work/` live on an external drive via `--data_dir` / `-w`.

## Status

Scaffold validated end-to-end (`-profile test -stub-run`, all 15 stages). Stage implementations, the comparative analysis layer, HPC packaging, and a real-SRA test are in progress — see `docs/`.
