# fungiforge — design decisions

- **Two layers, one contract.** Layer 1 (Nextflow, per-isolate) emits `master_fungi.tsv`
  + per-feature long tables; Layer 2 (R objective scripts) consumes only those. Mirrors
  the bacterial pipeline so the two studies are directly comparable and Layer 2 can be
  re-run independently.
- **ONT-first, hybrid-optional.** All isolates get Illumina eventually, but ONT lands
  first. Medaka always; Polypolish/POLCA when Illumina is present. Antifungal-resistance
  point-mutation calls are flagged **provisional on ONT-only assemblies** because
  homopolymer indels mimic the very substitutions/frameshifts we call.
- **Resistance is bespoke.** No ResFinder-for-fungi exists. A curated panel + FungAMR/MARDy
  references drive substitution/LoF/GOF calling, plus a dedicated structural detector for
  the *A. fumigatus* cyp51A promoter tandem repeats (TR34/TR46) — invisible to any SNP caller.
- **Identification without GTDB.** ITS/LSU (UNITE) + genome ANI (sourmash/skani) + MLST,
  combined by GCPSR multi-locus concordance. Species delimitation is reported with confidence,
  never as a bare label.
- **Containers per stage, base image for the rest.** Heavy stages use their upstream images
  (Funannotate, antiSMASH, dfam/tetools, Flye, Medaka, BUSCO); the light stages share
  `aduton1000/fungiforge`. Native linux/amd64 on HPC; emulated on the Apple-Silicon Mac.
- **Databases fetched up front, never mid-run** (`bin/fetch_references.sh`, aria2c-resumable),
  staged on an external drive via `--data_dir` because they total ~150–250 GB.
- **Organism-agnostic.** When the real fungal reads arrive, only inputs/params change — no code.
