# aduton1000/antismash-ff:8.0.0-r3 — antiSMASH 8 standalone, container-hardened for Nextflow.
#  (1) + procps: Nextflow's task wrapper needs `ps` to launch a containerised task; the
#      stock antismash/standalone:8.0.0 image omits it, so BGC tasks fail before running.
#  (2) world-writable module data dirs: antiSMASH lazily writes brawn alignment caches
#      (e.g. nrps_pks/data/nterm.fasta.brawn_cache) into its own install tree at runtime,
#      which fails as a non-root host UID (Nextflow runs -u $(id -u)). Making them writable
#      lets the cache be created on first use.
#  (3) pre-built brawn alignment caches at the path the NRPS/PKS order-finder actually reads.
#      antiSMASH 8.0.0 caches nterm/cterm alignments under nrps_pks/data/terminals/ but reads
#      them from nrps_pks/data/, so the first PKS cluster tries to write into the package dir —
#      fine on a writable Docker filesystem, fatal (read-only) under Apptainer/Singularity.
# Build: docker build --platform linux/amd64 -t aduton1000/antismash-ff:8.0.0-r3 -f env/antismash-ff.Dockerfile env
FROM antismash/standalone:8.0.0
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
RUN chmod -R a+rwX /usr/local/lib/python3.11/dist-packages/antismash || true
# (3) warm the brawn caches where orderfinder.perform_docking_domain_analysis looks for them
RUN python3 -c "import os; from antismash.common import brawn; import antismash.modules.nrps_pks.orderfinder as o; \
    d=os.path.join(os.path.dirname(o.__file__),'data'); \
    [brawn.get_cached_alignment(p, d) for p in (o.N_TERMINAL_PATH, o.C_TERMINAL_PATH)]; \
    print(sorted(f for f in os.listdir(d) if f.endswith('.brawn_cache')))" \
    && chmod -R a+rwX /usr/local/lib/python3.11/dist-packages/antismash || true
