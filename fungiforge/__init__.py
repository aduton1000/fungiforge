"""FungiForge — reproducible fungal genomics from ONT (+ optional Illumina hybrid).

Layer 1 is a Nextflow DSL2 pipeline (see main.nf); this Python package is the thin
CLI front door plus the bespoke per-isolate callers invoked by the pipeline's bin/
scripts (antifungal-resistance panel, novelty, report). Kept stdlib-light so the CLI
has no heavy dependencies — the real runtime tools live in the conda envs / containers.
"""
__version__ = "0.1.0"
