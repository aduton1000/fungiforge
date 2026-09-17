# fungiforge base image — the light tools for stages without a dedicated public image
# (identify / resistance / mge / extras-lite / report). Heavy stages use their own
# images (staphb/flye, ontresearch|staphb/medaka, dfam/tetools, nextgenusfs/funannotate,
# antismash/standalone, ezlabgva/busco). Built for linux/amd64:
#   docker build -f env/Dockerfile -t aduton1000/fungiforge:<version> .
# nQuire (ploidy from read allele frequencies; W2.6) is not packaged: built here from source with its
# bundled htslib (static, zlib only) and copied into the runtime image as a single binary.
FROM debian:bookworm-slim AS nquire
RUN apt-get update -q && apt-get install -y -q --no-install-recommends build-essential zlib1g-dev git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
ARG NQUIRE_COMMIT=master
RUN git clone --recursive -q https://github.com/clwgg/nQuire /opt/nQuire && cd /opt/nQuire && git checkout -q ${NQUIRE_COMMIT} \
    && make submodules && make && ./nQuire 2>&1 | head -1

FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.source="https://github.com/aduton1000/fungiforge"
USER root
COPY --from=nquire /opt/nQuire/nQuire /usr/local/bin/nQuire
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 MPLCONFIGDIR=/tmp/matplotlib
# Install from the explicit lock (exact package URLs + md5) so rebuilds reproduce the validated
# environment byte for byte; env/base.yml is the human-readable source of the same pins.
# Bootstrapping a NEW lock (after changing a pin): build once from the yml, validate, then
# export the lock from that image with bin/lock_env.sh:
#   docker build --build-arg ENV_SPEC=env/base.yml -f env/Dockerfile -t aduton1000/fungiforge:<ver> .
ARG ENV_SPEC=env/base.linux-64.lock
# (copied into a directory so the file keeps its extension: micromamba picks the parser from it)
COPY ${ENV_SPEC} /tmp/spec/
RUN spec=$(ls /tmp/spec/* | head -1) && micromamba install -y -n base -f "$spec" && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
# install the fungiforge python package (bespoke callers)
COPY pyproject.toml README.md /opt/fungiforge/
COPY fungiforge /opt/fungiforge/fungiforge
COPY bin /opt/fungiforge/bin
RUN pip install --no-cache-dir -e /opt/fungiforge && chmod +x /opt/fungiforge/bin/*.py 2>/dev/null || true
WORKDIR /work
