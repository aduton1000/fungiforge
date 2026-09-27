# fungiforge SignalP 6 image (W2.6, stage 13 secretome) — built ON SITE from the academic-licence
# package (https://services.healthtech.dtu.dk/services/SignalP-6.0/), never pushed to a public
# registry. Layers SignalP 6 (and CPU torch) onto the fungiforge base image, so every other extras
# tool is present too; a site config points the `extras` label at it.
#
#   docker build --platform linux/amd64 -f env/signalp6.Dockerfile \
#     --build-arg BASE=aduton1000/fungiforge:<version> \
#     --build-arg SIGNALP_TGZ=signalp-6.0i.fast.tar.gz \
#     -t aduton1000/fungiforge-signalp6:<version> <directory holding the tarball>
#
# (bin/hpc_install.sh --signalp <tarball> does this and converts the result to a .sif.)
# The "fast" package is the intended one; "slow_sequential" works the same way (~10x slower).
ARG BASE=aduton1000/fungiforge:0.2.0
FROM ${BASE}
ARG SIGNALP_TGZ=signalp-6.0i.fast.tar.gz
COPY ${SIGNALP_TGZ} /tmp/signalp6.tar.gz
RUN set -e; mkdir -p /tmp/sp6 && tar xzf /tmp/signalp6.tar.gz -C /tmp/sp6 \
    && PKG=$(find /tmp/sp6 -maxdepth 2 -type d -name 'signalp-6-package' | head -1) \
    && pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu torch \
    && pip install --no-cache-dir "$PKG" \
    && SPDIR=$(python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))") \
    && if [ -d "$PKG/models" ]; then mkdir -p "$SPDIR/model_weights" && cp -r "$PKG"/models/* "$SPDIR/model_weights/"; fi \
    && rm -rf /tmp/sp6 /tmp/signalp6.tar.gz \
    && signalp6 --help >/dev/null
