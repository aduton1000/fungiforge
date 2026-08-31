# fungiforge base image — the light tools for stages without a dedicated public image
# (identify / resistance / mge / extras-lite / report). Heavy stages use their own
# images (staphb/flye, ontresearch|staphb/medaka, dfam/tetools, nextgenusfs/funannotate,
# antismash/standalone, ezlabgva/busco). Build native linux/amd64 (no emulation):
#   docker build -f env/Dockerfile -t aduton1000/fungiforge:0.1.0 .
FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.source="https://github.com/aduton1000/fungiforge"
USER root
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 MPLCONFIGDIR=/tmp/matplotlib
COPY env/base.yml /tmp/base.yml
RUN micromamba install -y -n base -f /tmp/base.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
# install the fungiforge python package (bespoke callers)
COPY pyproject.toml README.md /opt/fungiforge/
COPY fungiforge /opt/fungiforge/fungiforge
COPY bin /opt/fungiforge/bin
RUN pip install --no-cache-dir -e /opt/fungiforge && chmod +x /opt/fungiforge/bin/*.py 2>/dev/null || true
WORKDIR /work
