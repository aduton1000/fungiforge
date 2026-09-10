#!/usr/bin/env bash
# ==============================================================================
# fungiforge — site-wide HPC installer (idempotent; rerun to update)
# ------------------------------------------------------------------------------
# Installs a central, read-only-for-users fungiforge under <prefix> in the same
# layout as callforge/captureforge on the same cluster, and wires the launchers
# onto PATH. Designed for a SLURM + Apptainer cluster WITHOUT a module system.
#
#   <prefix>/repo             git checkout (this repo; --ref selects tag/branch)
#   <prefix>/images/*.sif     the two fungiforge-built images (docker build -> SIF)
#   <prefix>/images/cache     shared Nextflow image cache (public images pre-pulled)
#   <prefix>/cli-env          conda env with the `fungiforge` Python CLI
#   <prefix>/site.config      nextflow -c overrides (caps, .sif paths)
#   <prefix>/fungiforge-env.sh   site env, sourced by the wrappers + /etc/profile.d
#   <prefix>/bin/{fungiforge,fungiforge-run,fungiforge-preflight,fungiforge-fetch-refs}
#   <db>                      reference DBs (filled by fungiforge-fetch-refs, hours)
#
# Usage (as a user who can write <prefix> and <db>, or with sudo for the mkdir step):
#   bash hpc_install.sh --prefix /hpc/opt/fungiforge --db /hpc/data/fungiforge --group hpcusers
# Options:
#   --prefix DIR      install root                         (default /hpc/opt/fungiforge)
#   --db DIR          reference database root              (default /hpc/data/fungiforge)
#   --group NAME      group that owns the install (users need read+exec) (default: primary group)
#   --repo URL        git remote                           (default github aduton1000/fungiforge)
#   --ref REF         tag/branch/commit to deploy          (default main)
#   --partition NAME  default SLURM partition for runs     (default: cluster default)
#   --bind PATHS      extra container bind roots           (default /hpc)
#   --skip-images     don't build the two .sif images / pre-pull public images
#   --rebuild-images  rebuild the .sif files even if present
#   --skip-cli        don't create the conda cli-env
#   --profile-d       write /etc/profile.d/fungiforge.sh via sudo (else prints the snippet)
#   --smoke           after install, run the stub-run DAG test through SLURM (~2-5 min)
#   --dry-run         print what would be done
# ==============================================================================
set -euo pipefail

PREFIX=/hpc/opt/fungiforge; DB=/hpc/data/fungiforge; GROUP=""; REPO=https://github.com/aduton1000/fungiforge.git
REF=main; PARTITION=""; BIND=/hpc; SKIP_IMAGES=0; REBUILD_IMAGES=0; SKIP_CLI=0; PROFILE_D=0; SMOKE=0; DRY=0
while [ $# -gt 0 ]; do case "$1" in
  --prefix) PREFIX="$2"; shift;;   --db) DB="$2"; shift;;   --group) GROUP="$2"; shift;;
  --repo) REPO="$2"; shift;;       --ref) REF="$2"; shift;; --partition) PARTITION="$2"; shift;;
  --bind) BIND="$2"; shift;;       --skip-images) SKIP_IMAGES=1;; --rebuild-images) REBUILD_IMAGES=1;;
  --skip-cli) SKIP_CLI=1;;         --profile-d) PROFILE_D=1;; --smoke) SMOKE=1;; --dry-run) DRY=1;;
  -h|--help) sed -n '2,36p' "$0"; exit 0;;
  *) echo "unknown option $1" >&2; exit 2;;
esac; shift; done
[ -n "$GROUP" ] || GROUP="$(id -gn)"

log(){ printf '\n\033[1;34m[install]\033[0m %s\n' "$*"; }
warn(){ printf '\033[1;33m[install] WARNING:\033[0m %s\n' "$*" >&2; }
die(){ printf '\033[1;31m[install] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }
run(){ if [ "$DRY" = 1 ]; then echo "  + $*"; else echo "  + $*"; "$@"; fi; }
has(){ command -v "$1" >/dev/null 2>&1; }

# ---- 0. preflight --------------------------------------------------------------
log "preflight"
for t in git nextflow java; do has "$t" || die "$t not on PATH"; done
RT=""; has apptainer && RT=apptainer; [ -z "$RT" ] && has singularity && RT=singularity
[ -n "$RT" ] || die "apptainer/singularity not on PATH"
NXF_V="$(nextflow -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
[ -n "$NXF_V" ] || die "nextflow on PATH ($(command -v nextflow)) did not report a version"
echo "  nextflow $NXF_V  |  $RT $($RT --version | awk '{print $NF}')  |  $(java -version 2>&1 | head -1)"
if [ "$SKIP_IMAGES" = 0 ]; then
  has docker && docker info >/dev/null 2>&1 || die "docker is needed to build the two fungiforge images (or pass --skip-images and build elsewhere)"
fi
if [ "$SKIP_CLI" = 0 ]; then
  CONDA=""; for c in mamba micromamba conda; do has "$c" && { CONDA="$c"; break; }; done
  [ -n "$CONDA" ] || die "mamba/conda needed for the cli-env (or --skip-cli)"
fi
for d in "$PREFIX" "$DB"; do
  if [ -d "$d" ]; then [ -w "$d" ] || die "$d exists but is not writable by $(id -un). Fix with:  sudo chown $(id -un):$GROUP '$d' && sudo chmod 2775 '$d'"
  else
    parent="$(dirname "$d")"
    [ -w "$parent" ] || die "cannot create $d ($parent not writable). Run once:  sudo mkdir -p '$d' && sudo chown $(id -un):$GROUP '$d' && sudo chmod 2775 '$d'"
  fi
done
case "$PREFIX$DB" in *' '*) die "install/db paths must not contain spaces (funannotate/antiSMASH break)";; esac

# ---- 1. layout + repo ---------------------------------------------------------
log "layout under $PREFIX (group $GROUP)"
run mkdir -p "$PREFIX/images/cache" "$PREFIX/bin" "$DB"
if [ -d "$PREFIX/repo/.git" ]; then
  log "updating repo -> $REF"
  run git -C "$PREFIX/repo" fetch --tags origin
  run git -C "$PREFIX/repo" checkout -q "$REF"
  git -C "$PREFIX/repo" symbolic-ref -q HEAD >/dev/null 2>&1 && run git -C "$PREFIX/repo" pull -q --ff-only || true
else
  log "cloning $REPO -> $PREFIX/repo"
  run git clone -q "$REPO" "$PREFIX/repo"
  run git -C "$PREFIX/repo" checkout -q "$REF"
fi
REPO_DIR="$PREFIX/repo"
VERSION="$(grep -oE "version *= *'[^']+'" "$REPO_DIR/nextflow.config" 2>/dev/null | head -1 | sed -E "s/.*'([^']+)'/\1/" || true)"
[ -n "$VERSION" ] || { [ "$DRY" = 1 ] && VERSION=0.1.0 || die "cannot read manifest.version from $REPO_DIR/nextflow.config"; }
COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "  fungiforge $VERSION @ $COMMIT"

# ---- 2. images ----------------------------------------------------------------
FF_TAG="aduton1000/fungiforge:$VERSION";      FF_SIF="$PREFIX/images/fungiforge-$VERSION.sif"
AS_TAG="aduton1000/antismash-ff:8.0.0-r2";    AS_SIF="$PREFIX/images/antismash-ff-8.0.0-r2.sif"
if [ "$SKIP_IMAGES" = 0 ]; then
  export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$PREFIX/images/.tmp}"; mkdir -p "$APPTAINER_TMPDIR"
  build_sif(){ # build_sif <tag> <dockerfile> <context> <sif>
    local tag="$1" df="$2" ctx="$3" sif="$4"
    if [ -s "$sif" ] && [ "$REBUILD_IMAGES" = 0 ]; then echo "  have $sif — skip (use --rebuild-images)"; return; fi
    log "docker build $tag  (native linux/amd64)"
    run docker build --platform linux/amd64 -t "$tag" -f "$df" "$ctx"
    log "convert -> $sif"
    # build to a temp name and rename only on success: an interrupted build (ssh drop,
    # Ctrl-C) must never leave a truncated .sif that a rerun would treat as finished.
    run rm -f "$sif" "$sif.part"
    run "$RT" build "$sif.part" "docker-daemon://$tag"
    run mv "$sif.part" "$sif"
  }
  build_sif "$FF_TAG" "$REPO_DIR/env/Dockerfile"              "$REPO_DIR"     "$FF_SIF"
  build_sif "$AS_TAG" "$REPO_DIR/env/antismash-ff.Dockerfile" "$REPO_DIR/env" "$AS_SIF"
  log "sanity: tools inside the images"
  run "$RT" exec "$FF_SIF" bash -c 'ps --version | head -1; fungiforge version; ITSx -h 2>&1 | head -1; sourmash --version'
  run "$RT" exec "$AS_SIF" bash -c 'ps --version | head -1; antismash --version'
  log "pre-pull public images into the shared cache (funannotate ≈15 GB — be patient)"
  export FUNGIFORGE_DB="$DB" NXF_APPTAINER_CACHEDIR="$PREFIX/images/cache" NXF_SINGULARITY_CACHEDIR="$PREFIX/images/cache"
  if [ "$DRY" = 0 ]; then
    rm -f "$DB/containers/.done"
    bash "$REPO_DIR/bin/fetch_references.sh" images || warn "some public images failed to pull — see $DB/logs/images.log; Nextflow will retry at run time"
  fi
fi

# ---- 3. cli-env ---------------------------------------------------------------
if [ "$SKIP_CLI" = 0 ]; then
  if [ -x "$PREFIX/cli-env/bin/python" ]; then
    log "cli-env exists — reinstalling the fungiforge package"
  else
    log "creating $PREFIX/cli-env ($CONDA)"
    export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$HOME/.conda/pkgs}"
    run "$CONDA" create -y -q -p "$PREFIX/cli-env" -c conda-forge python=3.11 pip
  fi
  run "$PREFIX/cli-env/bin/pip" install -q --no-cache-dir "$REPO_DIR"
  run "$PREFIX/cli-env/bin/fungiforge" version
fi

# ---- 4. site.config + env + wrappers ------------------------------------------
log "site.config / fungiforge-env.sh / wrappers"
if [ "$DRY" = 0 ]; then
  for need in share/site.config.example share/fungiforge-env.sh.example share/bin/fungiforge-run bin/fetch_references.sh bin/preflight_qc.sh; do
    [ -f "$REPO_DIR/$need" ] || die "$REPO_DIR/$need missing — the deployed ref '$REF' predates the HPC kit; push the latest main (or pass --ref) and rerun"
  done
  if [ ! -f "$PREFIX/site.config" ]; then
    sed -e "s#/hpc/opt/fungiforge#$PREFIX#g" \
        -e "s#fungiforge-0.1.0.sif#$(basename "$FF_SIF")#" \
        "$REPO_DIR/share/site.config.example" > "$PREFIX/site.config"
    if [ -n "$PARTITION" ]; then   # (no sed -i: differs between GNU and BSD)
      sed "s#// slurm_partition = 'global'.*#slurm_partition = '$PARTITION'#" "$PREFIX/site.config" > "$PREFIX/site.config.tmp" && mv "$PREFIX/site.config.tmp" "$PREFIX/site.config"
    fi
  else echo "  keeping existing $PREFIX/site.config"; fi
  if [ ! -f "$PREFIX/fungiforge-env.sh" ]; then
    sed -e "s#^export FUNGIFORGE_ROOT=.*#export FUNGIFORGE_ROOT=\"$PREFIX\"#" \
        -e "s#^export FUNGIFORGE_DB=.*#export FUNGIFORGE_DB=\"$DB\"#" \
        -e "s#^export FUNGIFORGE_PARTITION=.*#export FUNGIFORGE_PARTITION=\"$PARTITION\"#" \
        -e "s#^export FUNGIFORGE_BIND=.*#export FUNGIFORGE_BIND=\"$BIND\"#" \
        "$REPO_DIR/share/fungiforge-env.sh.example" > "$PREFIX/fungiforge-env.sh"
  else echo "  keeping existing $PREFIX/fungiforge-env.sh"; fi

  cat > "$PREFIX/bin/fungiforge-run" <<EOF
#!/usr/bin/env bash
source "$PREFIX/fungiforge-env.sh"
exec bash "\$FUNGIFORGE_HOME/share/bin/fungiforge-run" "\$@"
EOF
  cat > "$PREFIX/bin/fungiforge" <<EOF
#!/usr/bin/env bash
source "$PREFIX/fungiforge-env.sh"
exec "$PREFIX/cli-env/bin/fungiforge" "\$@"
EOF
  cat > "$PREFIX/bin/fungiforge-preflight" <<EOF
#!/usr/bin/env bash
# pre-run ONT triage (docs/preflight_qc.md); runs inside the fungiforge image
source "$PREFIX/fungiforge-env.sh"
export FUNGIFORGE_SIF="$FF_SIF"
exec bash "\$FUNGIFORGE_HOME/bin/preflight_qc.sh" "\$@"
EOF
  cat > "$PREFIX/bin/fungiforge-fetch-refs" <<EOF
#!/usr/bin/env bash
# one-time reference-database download into \$FUNGIFORGE_DB (hours; resumable)
source "$PREFIX/fungiforge-env.sh"
export NXF_APPTAINER_CACHEDIR="\$FUNGIFORGE_IMAGE_CACHE" NXF_SINGULARITY_CACHEDIR="\$FUNGIFORGE_IMAGE_CACHE"
exec bash "\$FUNGIFORGE_HOME/bin/fetch_references.sh" "\$@"
EOF
  chmod 755 "$PREFIX"/bin/* "$REPO_DIR"/bin/*.sh "$REPO_DIR"/share/bin/* 2>/dev/null || true
fi

# ---- 5. permissions -----------------------------------------------------------
log "permissions: group $GROUP read/exec; images/cache group-writable (setgid)"
if [ "$DRY" = 0 ]; then
  chgrp -R "$GROUP" "$PREFIX" 2>/dev/null || warn "chgrp failed on some files (not owner) — run: sudo chgrp -R $GROUP $PREFIX"
  chmod -R u+rwX,g+rX,o+rX "$PREFIX" 2>/dev/null || true
  chmod 2775 "$PREFIX/images/cache" 2>/dev/null || true
  chmod -R g+rwX "$PREFIX/images/cache" 2>/dev/null || true
  [ -d "$DB" ] && { chgrp "$GROUP" "$DB" 2>/dev/null || true; chmod 2775 "$DB" 2>/dev/null || true; }
fi

# ---- 6. PATH hook -------------------------------------------------------------
SNIPPET="# fungiforge — site-wide launchers (installed by bin/hpc_install.sh)
[ -r $PREFIX/fungiforge-env.sh ] && source $PREFIX/fungiforge-env.sh"
if [ "$PROFILE_D" = 1 ] && [ "$DRY" = 0 ]; then
  log "writing /etc/profile.d/fungiforge.sh (sudo)"
  echo "$SNIPPET" | sudo tee /etc/profile.d/fungiforge.sh >/dev/null && sudo chmod 644 /etc/profile.d/fungiforge.sh
else
  log "PATH hook — add this to /etc/profile.d/fungiforge.sh (all users) or ~/.bashrc:"
  echo "$SNIPPET"
fi

# ---- 7. smoke test ------------------------------------------------------------
if [ "$SMOKE" = 1 ] && [ "$DRY" = 0 ]; then
  log "smoke test: 15-stage DAG stub-run through SLURM"
  SM="$HOME/fungiforge_smoke_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$SM"; cd "$SM"
  source "$PREFIX/fungiforge-env.sh"
  bash "$REPO_DIR/share/bin/fungiforge-run" -profile "test,hpc_slurm,$RT" -stub-run --outdir "$SM/results" \
    && log "SMOKE OK — $SM/results/pipeline_info/execution_report.html" \
    || die "smoke test failed — see $SM/.nextflow.log"
fi

log "done."
cat <<EOF
  install root : $PREFIX   (fungiforge $VERSION @ $COMMIT)
  databases    : $DB   -> run  fungiforge-fetch-refs   (one-time, hours; resumable)
  images       : $FF_SIF
                 $AS_SIF
                 cache: $PREFIX/images/cache
  next         : source $PREFIX/fungiforge-env.sh; fungiforge-run --help
EOF
