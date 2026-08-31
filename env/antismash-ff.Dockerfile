# aduton1000/antismash-ff:8.0.0 — antiSMASH 8 standalone, container-hardened for Nextflow.
#  (1) + procps: Nextflow's task wrapper needs `ps` to launch a containerised task; the
#      stock antismash/standalone:8.0.0 image omits it, so BGC tasks fail before running.
#  (2) world-writable module data dirs: antiSMASH lazily writes brawn alignment caches
#      (e.g. nrps_pks/data/nterm.fasta.brawn_cache) into its own install tree at runtime,
#      which fails as a non-root host UID (Nextflow runs -u $(id -u)). Making them writable
#      lets the cache be created on first use.
# Build: docker build --platform linux/amd64 -t aduton1000/antismash-ff:8.0.0-r2 -f env/antismash-ff.Dockerfile env
FROM antismash/standalone:8.0.0
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
RUN chmod -R a+rwX /usr/local/lib/python3.11/dist-packages/antismash || true
