# fungiforge — design decisions

- **Two layers, one contract.** Layer 1 (Nextflow, per-isolate) emits `master_fungi.tsv`
  + per-feature long tables; Layer 2 (R objective scripts) consumes only those. Mirrors
  the bacterial pipeline so the two studies are directly comparable and Layer 2 can be
  re-run independently.
- **Platform-flexible, one contract.** Each isolate is routed by its reads
  (`meta.assembly_mode`): **long-read** (ONT present) uses Flye + Medaka, with Polypolish
  adding Illumina when it's there (hybrid); **short-read** (Illumina only, no ONT) uses
  SPAdes. All three modes converge on the same `*.polished.fasta` → downstream stages, so
  nothing after Stage 03b needs to know how the assembly was made. Antifungal-resistance
  point-mutation calls are flagged **provisional on ONT-only assemblies** (homopolymer
  indels mimic the very substitutions/frameshifts we call) but **high-confidence on hybrid
  and Illumina-only assemblies**, which have no homopolymer-indel problem. Illumina-only
  isolates skip Medaka/Polypolish (nothing to hybrid-polish) and report `mode=illumina_only`.
- **Resistance is bespoke.** No ResFinder-for-fungi exists. A curated panel + FungAMR/MARDy
  references drive substitution/LoF/GOF calling, plus a dedicated structural detector for
  the *A. fumigatus* cyp51A promoter tandem repeats (TR34/TR46) — invisible to any SNP caller.
- **Identification without GTDB.** ITS/LSU (UNITE) + genome ANI (sourmash/skani) + MLST,
  combined by GCPSR multi-locus concordance. Species delimitation is reported with confidence,
  never as a bare label.
- **Containers per stage, base image for the rest.** Heavy stages use their upstream images
  (Funannotate, antiSMASH, dfam/tetools, Flye, Medaka, BUSCO); the light stages share
  `aduton1000/fungiforge`. Native linux/amd64 on x86_64 Linux; run under emulation on arm64 hosts.
- **Databases fetched up front, never mid-run** (`bin/fetch_references.sh`, aria2c-resumable),
  staged on an external drive via `--data_dir` because they total ~150–250 GB.
- **Organism-agnostic.** When the real fungal reads arrive, only inputs/params change — no code.
- **Stay on funannotate v1 for now; funannotate2 is not a drop-in (decided 2026-09-20, W6.2).**
  Evaluated against the released sources rather than the changelog. funannotate2 (v26.6.21,
  2026-06-22) would *remove* four capabilities this pipeline depends on. It has **no Trinity/PASA
  RNA-seq assembly path** — RNA enters `predict` only as an aligned BAM — so it cannot deliver
  W6.2's own goal. It has **no `annotate_results/*.annotations.txt`**: functional results move to
  three-column TSVs under `annotate_misc/`, which breaks `annotate_stats.py` and the four coverage
  columns. It has **no `--eggnog` / `--iprscan`**, so stages 07b and 07c would need reworking into
  its `-a/--annotations` format. And it drops `update` entirely, so there is no UTR/model
  correction step. Packaging is also not ready: Docker Hub carries no tag for the current release
  and the BioContainers build is Python-only, without Augustus or GeneMark. Neither project states
  that v1 is deprecated, and v1's master branch is still developed. **Revisit when** funannotate2
  gains a transcript-assembly path or a documented equivalent of the wide annotations table, and a
  released, predictor-complete image exists. Until then v1.8.17 stays pinned by digest.
- **`meta` is a cache key, not a scratchpad (2026-09-20).** Every Nextflow task's signature includes the
  per-sample `meta` map, so adding one key invalidates every cached task of every existing run. Adding
  `has_rna` for W6.2 cost a completed CEA10 validation its entire cache and would have rerun ten hours of
  work. Keep `meta` to `id`, `assembly_mode`, `compartment`, `facility`, `season`; express anything a single
  stage needs through the files staged into that stage, which is where the pipeline already answers it.
