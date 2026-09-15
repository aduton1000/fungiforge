# fungiforge base image — the light tools for stages without a dedicated public image
# (identify / resistance / mge / extras-lite / report). Heavy stages use their own
# images (staphb/flye, ontresearch|staphb/medaka, dfam/tetools, nextgenusfs/funannotate,
# antismash/standalone, ezlabgva/busco). Built for linux/amd64:
#   docker build -f env/Dockerfile -t aduton1000/fungiforge:<version> .
FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.source="https://github.com/aduton1000/fungiforge"
USER root
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
