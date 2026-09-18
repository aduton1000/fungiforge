#!/usr/bin/env bash
# ==============================================================================
# fungiforge — fetch reference databases + container images
# ------------------------------------------------------------------------------
# Downloads every database fungiforge needs into a single data dir (the forge
# convention: fetch + checksum UP FRONT, never mid-run). Idempotent and
# resumable — each DB is guarded by a `.done` sentinel, and one DB failing does
# NOT abort the others. Designed to run for many hours on an external SSD.
#
# Usage:
#   export FUNGIFORGE_DB="/path/to/fungiforge_db"
#   bin/fetch_references.sh                 # all steps
#   bin/fetch_references.sh images antismash funannotate   # selected steps
#
# Steps: images antismash funannotate eggnog busco unite kraken2 refseq_fungi
#        fungamr rvdb benchmarks markers mlst dbcan phibase effectorp genomad genomes interproscan
#
# NB several DBs are downloaded THROUGH their tool container (funannotate setup,
# eggnog, antismash) so the relevant image is pulled first. Version/URL-sensitive
# DBs (UNITE, Kraken2, sourmash-fungi) resolve the CURRENT release at run time;
# if a URL 404s the step logs and skips rather than guessing.
# ==============================================================================
set -uo pipefail

DB="${FUNGIFORGE_DB:?set FUNGIFORGE_DB to the target data dir, e.g. /data/fungiforge_db}"
mkdir -p "$DB"/{containers,funannotate,eggnog,antismash,busco,unite,kraken2,refseq_fungi,fungamr,rvdb,markers,mlst,interproscan,dbcan,phibase,effectorp,logs}
LOGDIR="$DB/logs"; MANIFEST="$DB/MANIFEST.tsv"

# One fetch at a time per data dir. Two concurrent runs resume the same partial file and
# truncate each other's work — a download that grows, shrinks and never finishes.
if command -v flock >/dev/null 2>&1; then
  exec 9>"$DB/.fetch.lock"
  if ! flock -n 9; then
    echo "another fetch_references run is using $DB (lock: $DB/.fetch.lock)." >&2
    # Name the holder. Killing the script alone can leave a download running that still holds the
    # lock, and "pkill -f fetch_references" then looks like it did nothing.
    holder="$( { fuser "$DB/.fetch.lock" 2>/dev/null || lsof -t "$DB/.fetch.lock" 2>/dev/null; } | tr -s ' ' '\n' | grep -E '^[0-9]+$' | sort -u | tr '\n' ' ')"
    if [ -n "${holder// /}" ]; then
      echo "held by:" >&2; ps -o pid=,etime=,args= -p ${holder} 2>/dev/null | sed 's/^/  /' >&2
      echo "wait for it, or stop it:  kill ${holder}" >&2
    else
      echo "wait for it, or stop it first:  pkill -f fetch_references" >&2
    fi
    exit 1
  fi
fi
[ -f "$MANIFEST" ] || echo -e "database\tdetail\tstatus\ttimestamp" > "$MANIFEST"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/fetch.log" ; }

# The downloader must not inherit the lock file descriptor (9<&- below): a download that outlives a
# killed parent would otherwise keep the data dir locked, and the next run reports a fetch in
# progress that no longer exists.
# Robust, resumable, multi-connection download. aria2c preferred (16 parallel
# connections, --continue resumes a broken transfer from where it stopped, infinite
# retries); curl -C - is the fallback. Never restarts a partial file from scratch.
dl(){  # dl <url> <output-filename>   (run inside the target dir)
  local url="$1" out="$2"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c --continue=true -x16 -s16 -k1M --max-tries=0 --retry-wait=10 \
           --file-allocation=none --auto-file-renaming=false --console-log-level=warn \
           -o "$out" "$url" 9<&-
  elif command -v wget >/dev/null 2>&1; then
    # Preferred over curl for the multi-GB databases. curl's --retry rewinds to byte 0 when the
    # server drops a transfer that -C - had positioned, so a flaky link makes the file grow and
    # shrink forever; wget -c re-issues a Range request from the bytes already on disk on every
    # retry. -nv keeps one line per event in the log instead of a rewriting progress bar.
    wget -c -nv --tries=0 --timeout=60 --read-timeout=300 --waitretry=30 -O "$out" "$url" 9<&-
  else
    # -sS: no progress meter. Everything here is appended to a log file, where curl's meter
    # rewrites one line thousands of times and buries the real messages; errors still print.
    # Watch progress by file size instead:  watch -n 20 'ls -lh <db>/<dir>'
    curl -fL -C - --retry 999 --retry-delay 10 --retry-all-errors -sS -o "$out" "$url" 9<&-
  fi
}

# A resumed multi-GB download that picked up from a corrupt offset still decompresses to a short,
# useless database. Check the archive before trusting it; a failure here deletes the file so the
# next run starts clean rather than resuming onto the same damage.
verify_archive(){  # verify_archive <file>   (run inside the target dir)
  local f="$1"
  case "$f" in
    *.tar.gz|*.tgz) tar tzf "$f" >/dev/null 2>&1 ;;
    *.gz)     gzip -t "$f" 2>/dev/null ;;
    *)        return 0 ;;
  esac || { echo "[fetch] $f is corrupt (failed its integrity check) — removing it" >&2; rm -f "$f"; return 1; }
}


# Print the growing file's size every 60 s while a download runs in the background, so a long
# fetch shows progress in the log without a progress meter. Usage: watch_size <file> & ; WPID=$!
watch_size(){
  local f="$1" last=0 now
  while sleep 60; do
    [ -e "$f" ] || continue
    now=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f" 2>/dev/null || echo 0)
    [ "$now" = "$last" ] && continue
    log "  $(basename "$f"): $(numfmt --to=iec "$now" 2>/dev/null || echo "$now")"
    last=$now
  done
}

# ---- run a command inside a container image, whatever runtime the host has (L24) -------------
# crun <image> <host_dir>:<container_dir> [<host_dir>:<container_dir> ...] -- <command...>
# Apptainer/Singularity first (the cluster path: no daemon, no root), then Docker. Apptainer
# ignores an image ENTRYPOINT with `exec`, which is what both callers below need; under Docker the
# entrypoint is cleared explicitly. Images are cached under $DB/containers for the apptainer path.
crun(){
  local img="$1"; shift
  local binds=() mounts=()
  while [ "${1:-}" != "--" ] && [ $# -gt 0 ]; do binds+=("$1"); shift; done
  [ "${1:-}" = "--" ] && shift
  local rt=""
  command -v apptainer >/dev/null 2>&1 && rt=apptainer
  [ -z "$rt" ] && command -v singularity >/dev/null 2>&1 && rt=singularity
  if [ -n "$rt" ]; then
    local cache="${NXF_APPTAINER_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR:-$DB/containers}}"
    mkdir -p "$cache"
    local sif
    sif="$cache/$(echo "$img" | sed -E 's#[/:]#-#g').img"
    if [ ! -s "$sif" ]; then
      log "  $rt pull $img"
      rm -f "$sif.part"
      $rt pull --name "$sif.part" "docker://$img" >>"$LOGDIR/containers.log" 2>&1 || return 1
      mv "$sif.part" "$sif"
    fi
    for b in "${binds[@]}"; do mounts+=(-B "$b"); done
    $rt exec "${mounts[@]}" "$sif" "$@"
  elif command -v docker >/dev/null 2>&1; then
    for b in "${binds[@]}"; do mounts+=(-v "$b"); done
    docker run --rm --platform linux/amd64 --entrypoint "" "${mounts[@]}" "$img" "$@"
  else
    echo "no apptainer/singularity/docker on PATH" >&2; return 127
  fi
}

is_done(){ [ -f "$DB/$1/.done" ]; }
mark(){ date '+%F %T' > "$DB/$1/.done"; echo -e "$1\t$2\tOK\t$(date '+%F %T')" >> "$MANIFEST"; log "✓ $1 done ($2)"; }
fail(){ echo -e "$1\t$2\tFAILED\t$(date '+%F %T')" >> "$MANIFEST"; log "✗ $1 FAILED ($2) — continuing"; }

# ---- container images (linux/amd64; native on x86_64 Linux, emulated on arm64 hosts) ----
# Image list is derived from conf/base.config so it never drifts. The two
# aduton1000/* images are built locally (env/Dockerfile, env/antismash-ff.Dockerfile)
# and are skipped here. On a Docker host: `docker pull`. On an Apptainer/Singularity
# host: pre-pull into the shared Nextflow image cache ($NXF_APPTAINER_CACHEDIR /
# $NXF_SINGULARITY_CACHEDIR, default $DB/containers) using Nextflow's own file-naming
# (registry/repo:tag -> registry-repo-tag.img) so runs find them without pulling.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
public_images(){
  grep -oE "container *= *'[^']+'" "$REPO_DIR/conf/base.config" | sed -E "s/.*'([^']+)'/\1/" \
    | grep -v '^aduton1000/' | sort -u
}
step_images(){
  is_done containers && { log "images already pulled — skip"; return; }
  local ok=1 img
  local rt=""; command -v apptainer >/dev/null 2>&1 && rt=apptainer
  [ -z "$rt" ] && command -v singularity >/dev/null 2>&1 && rt=singularity
  if [ -n "$rt" ]; then
    local cache="${NXF_APPTAINER_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR:-$DB/containers}}"
    mkdir -p "$cache"
    for img in $(public_images); do
      local name; name="$(echo "$img" | sed -E 's#[/:]#-#g').img"
      if [ -s "$cache/$name" ]; then log "have $name — skip"; continue; fi
      log "$rt pull $img -> $cache/$name"
      # pull to a temp name, rename on success (an interrupted pull must not look finished)
      rm -f "$cache/$name.part"
      { $rt pull --name "$cache/$name.part" "docker://$img" >>"$LOGDIR/images.log" 2>&1 && mv "$cache/$name.part" "$cache/$name"; } \
        || { fail containers "$img"; ok=0; rm -f "$cache/$name.part"; }
    done
    [ "$ok" = 1 ] && mark containers "$(public_images | wc -l | tr -d ' ') images -> $cache"
  elif command -v docker >/dev/null 2>&1; then
    for img in $(public_images); do
      log "docker pull $img"
      docker pull --platform linux/amd64 "$img" >>"$LOGDIR/images.log" 2>&1 || { fail containers "$img"; ok=0; }
    done
    [ "$ok" = 1 ] && mark containers "$(public_images | wc -l | tr -d ' ') images (docker)"
  else
    fail containers "no apptainer/singularity/docker on PATH"
  fi
}

# ---- antiSMASH databases (~9 GB) via the antismash 8 image -------------------
step_antismash(){
  is_done antismash && { log "antismash db present — skip"; return; }
  log "downloading antiSMASH databases -> $DB/antismash"
  # The image entrypoint is `antismash`; the DB downloader is a separate console script, so the
  # entrypoint is bypassed (crun does that for both runtimes). It takes --database-dir.
  crun antismash/standalone:8.0.0 "$DB/antismash:/db" -- download-antismash-databases --database-dir /db >>"$LOGDIR/antismash.log" 2>&1 \
    && mark antismash "antismash8 db" || fail antismash "download-antismash-databases"
}

# ---- Funannotate database (~30-50 GB): Pfam, dbCAN, MEROPS, InterPro, BUSCO --
step_funannotate(){
  is_done funannotate && { log "funannotate db present — skip"; return; }
  log "funannotate setup -i all -> $DB/funannotate (large, hours)"
  FUNANNOTATE_DB=/data crun nextgenusfs/funannotate:latest "$DB/funannotate:/data" -- funannotate setup -i all -d /data >>"$LOGDIR/funannotate.log" 2>&1 \
    && mark funannotate "funannotate setup all" || fail funannotate "funannotate setup"
}

# ---- eggNOG database (~50 GB) -----------------------------------------------
# eggnog-mapper ships download_eggnog_data.py; run it from the funannotate image
# (bundles eggnog-mapper) or the eggnog-mapper image. -y = accept all downloads.
step_eggnog(){
  # eggNOG 5 data for eggNOG-mapper 2.1 (stage 07b): direct download, no container needed
  # (the funannotate image ships no emapper; stage 07b runs the eggnog-mapper image).
  # A file is only accepted when it reaches its expected size: a truncated leftover from an
  # interrupted or failed earlier attempt is non-empty, and treating that as done shipped a
  # 1.5 GB eggnog directory that emapper cannot use.
  is_done eggnog && { log "eggnog db present — skip"; return; }
  local base="${EGGNOG_URL:-http://eggnog5.embl.de/download/emapperdb-5.0.2}"
  local ok=1 f name min have
  log "downloading eggNOG 5.0.2 data (eggnog.db ~12 GB and eggnog_proteins.dmnd ~9 GB once decompressed) -> $DB/eggnog"
  for f in eggnog.db.gz:eggnog.db:8000000000 eggnog_proteins.dmnd.gz:eggnog_proteins.dmnd:4000000000 eggnog.taxa.tar.gz:eggnog.taxa.db:1000000; do
    name="${f#*:}"; min="${name#*:}"; name="${name%%:*}"; f="${f%%:*}"
    have=$(stat -c%s "$DB/eggnog/$name" 2>/dev/null || stat -f%z "$DB/eggnog/$name" 2>/dev/null || echo 0)
    if [ "$have" -ge "$min" ]; then log "have $name ($(numfmt --to=iec "$have" 2>/dev/null || echo "$have")) — skip"; continue; fi
    [ "$have" -gt 0 ] && log "  $name is $(numfmt --to=iec "$have" 2>/dev/null || echo "$have"), expected >= $(numfmt --to=iec "$min" 2>/dev/null || echo "$min") — re-downloading"
    ( cd "$DB/eggnog" && rm -f "$name" && dl "$base/$f" "$f" && verify_archive "$f" \
      && { case "$f" in *.tar.gz) tar xzf "$f" && rm -f "$f";; *.gz) gunzip -f "$f";; esac; } ) >>"$LOGDIR/eggnog.log" 2>&1 \
      || { fail eggnog "$f"; ok=0; }
  done
  local db_size dmnd_size
  db_size=$(stat -c%s "$DB/eggnog/eggnog.db" 2>/dev/null || stat -f%z "$DB/eggnog/eggnog.db" 2>/dev/null || echo 0)
  dmnd_size=$(stat -c%s "$DB/eggnog/eggnog_proteins.dmnd" 2>/dev/null || stat -f%z "$DB/eggnog/eggnog_proteins.dmnd" 2>/dev/null || echo 0)
  if [ "$ok" = 1 ] && [ "$db_size" -ge 8000000000 ] && [ "$dmnd_size" -ge 4000000000 ]; then
    mark eggnog "emapperdb-5.0.2 (eggnog.db $(numfmt --to=iec "$db_size" 2>/dev/null || echo "$db_size"), eggnog_proteins.dmnd $(numfmt --to=iec "$dmnd_size" 2>/dev/null || echo "$dmnd_size"), taxa)"
  else
    fail eggnog "incomplete after download (eggnog.db=$db_size dmnd=$dmnd_size); see logs/eggnog.log"
  fi
}

# ---- InterProScan data release (W2.5, stage 07c) ------------------------------
# The image (interpro/interproscan:<ver>, conf/base.config) carries the software only; the member
# databases come from the matching data tarball (~6.9 GB), extracted to
# $DB/interproscan/interproscan-<ver>/data and bound at /opt/interproscan/data by the profile
# (--interproscan_data). IPS_VERSION must match the image tag.
step_interproscan(){
  is_done interproscan && { log "interproscan data present — skip"; return; }
  local ver="${IPS_VERSION:-5.78-109.0}"
  local url="${IPS_DATA_URL:-https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$ver/alt/interproscan-data-$ver.tar.gz}"
  mkdir -p "$DB/interproscan"
  log "downloading InterProScan data $ver -> $DB/interproscan"
  ( cd "$DB/interproscan" && dl "$url" "interproscan-data-$ver.tar.gz" && dl "$url.md5" "interproscan-data-$ver.tar.gz.md5" \
    && md5sum -c "interproscan-data-$ver.tar.gz.md5" && tar xzf "interproscan-data-$ver.tar.gz" && rm -f "interproscan-data-$ver.tar.gz" ) >>"$LOGDIR/interproscan.log" 2>&1 \
    && [ -d "$DB/interproscan/interproscan-$ver/data" ] && mark interproscan "$ver -> interproscan-$ver/data" || fail interproscan "$url"
}

# ---- BUSCO / compleasm fungal lineages --------------------------------------
# compleasm is fastest; if not installed, fall back to BUSCO container download.
step_busco(){
  # fungi_odb10 is the lineage of the early QC gate (stage 05); the others are the species-aware
  # lineages stage 08b re-runs BUSCO with (fungiforge/resources/busco_lineages.tsv). Override the
  # set with BUSCO_LINEAGES="fungi_odb10 eurotiales_odb10 ..." (space-separated). Lineages are
  # downloaded straight from the BUSCO data server (no container needed) into
  # $DB/busco/lineages/<lineage>, the layout `busco --download_path $DB/busco --offline` and
  # compleasm `-L $DB/busco/lineages` read; file_versions.tsv is refreshed alongside.
  is_done busco && { log "busco lineages present — skip"; return; }
  local lineages="${BUSCO_LINEAGES:-fungi_odb10 eurotiales_odb10 saccharomycetes_odb10 hypocreales_odb10 tremellomycetes_odb10 mucorales_odb10 onygenales_odb10 sordariomycetes_odb10}"
  local base="${BUSCO_DATA_URL:-https://busco-data.ezlab.org/v5/data}"
  local ok=1 l ver
  mkdir -p "$DB/busco/lineages"
  log "fetching BUSCO lineages ($lineages) -> $DB/busco/lineages"
  ( cd "$DB/busco" && dl "$base/file_versions.tsv" file_versions.tsv ) >>"$LOGDIR/busco.log" 2>&1 \
    || { fail busco "file_versions.tsv from $base"; return; }
  for l in $lineages; do
    if [ -d "$DB/busco/lineages/$l" ] || [ -d "$DB/busco/$l" ] || [ -d "$DB/busco/busco_downloads/lineages/$l" ]; then log "have $l — skip"; continue; fi
    ver=$(awk -v L="$l" -F'\t' '$1==L {print $2; exit}' "$DB/busco/file_versions.tsv")
    if [ -z "$ver" ]; then fail busco "$l not in file_versions.tsv"; ok=0; continue; fi
    log "downloading $l ($ver)"
    ( cd "$DB/busco/lineages" && dl "$base/lineages/$l.$ver.tar.gz" "$l.tar.gz" && verify_archive "$l.tar.gz" && tar xzf "$l.tar.gz" && rm -f "$l.tar.gz" ) >>"$LOGDIR/busco.log" 2>&1 \
      && [ -d "$DB/busco/lineages/$l" ] || { fail busco "download/extract $l"; ok=0; }
  done
  [ "$ok" = 1 ] && mark busco "$lineages"
}

# ---- UNITE fungal ITS reference (version-sensitive URL) ----------------------
# The current general FASTA release URL changes per version. To refresh: find the
# latest DOI at https://unite.ut.ee/repository.php, then resolve the file link via
#   curl -s "https://api.plutof.ut.ee/v1/public/dois/?identifier=<DOI>" | python3 -c \
#     'import sys,json;print(json.load(sys.stdin)["data"][0]["attributes"]["media"][0]["url"])'
# Default below = UNITE general FASTA release for Fungi v10.0 (19.02.2025),
# DOI 10.15156/BIO/3301229.
step_unite(){
  is_done unite && { log "unite present — skip"; return; }
  local UNITE_URL="${UNITE_URL:-https://s3.hpc.ut.ee/plutof-public/original/9489f7bc-7cc1-4e0a-84dc-c732476b9acd.tgz}"
  if [ -z "$UNITE_URL" ]; then fail unite "UNITE_URL not set — resolve current release"; return; fi
  log "downloading UNITE -> $DB/unite"
  ( cd "$DB/unite" && dl "$UNITE_URL" unite.tgz && verify_archive unite.tgz && tar xzf unite.tgz ) >>"$LOGDIR/unite.log" 2>&1 \
    && mark unite "$UNITE_URL" || fail unite "$UNITE_URL"
}

# ---- Kraken2 contamination DB (PlusPF-8 by default; version-sensitive) -------
step_kraken2(){
  is_done kraken2 && { log "kraken2 present — skip"; return; }
  # Default = Kraken2 PlusPF-8GB, 2025-04-02 build (latest at genome-idx.s3 as of
  # 2026-07). Newer builds: list `curl -s "https://genome-idx.s3.amazonaws.com/?list-type=2&prefix=kraken/k2_pluspf_08gb"`.
  local K2_URL="${K2_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20250402.tar.gz}"
  if [ -z "$K2_URL" ]; then fail kraken2 "K2_URL not set — resolve current index at genome-idx.s3"; return; fi
  log "downloading Kraken2 DB -> $DB/kraken2"
  ( cd "$DB/kraken2" && dl "$K2_URL" k2.tgz && verify_archive k2.tgz && tar xzf k2.tgz ) >>"$LOGDIR/kraken2.log" 2>&1 \
    && mark kraken2 "$K2_URL" || fail kraken2 "$K2_URL"
}

# ---- RefSeq fungi genomes for sourmash gather / skani (version-sensitive) -----
# Prefer a prepared sourmash fungal database (OSF) if resolvable; else NCBI datasets.
step_refseq_fungi(){
  is_done refseq_fungi && { log "refseq_fungi present — skip"; return; }
  # Default = prepared sourmash GenBank fungi k=31 DB (genbank-2022.03, 10,286 genomes).
  local SM_URL="${SOURMASH_FUNGI_URL:-https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2022.03/genbank-2022.03-fungi-k31.zip}"
  if [ -n "$SM_URL" ]; then
    log "downloading prepared sourmash fungal DB -> $DB/refseq_fungi"
    ( cd "$DB/refseq_fungi" && dl "$SM_URL" fungi-k31.zip ) >>"$LOGDIR/refseq_fungi.log" 2>&1 \
      && mark refseq_fungi "$SM_URL" || fail refseq_fungi "$SM_URL"
  else
    fail refseq_fungi "SOURMASH_FUNGI_URL not set — resolve prepared DB or build from NCBI datasets"
  fi
}

# ---- FungAMR + MARDy antifungal-resistance catalogs (small) -------------------
step_fungamr(){
  is_done fungamr && { log "fungamr present — skip"; return; }
  log "fetching FungAMR antifungal-resistance catalogs -> $DB/fungamr"
  # FungAMR (Bédard/Landry et al., Nat Microbiol 2025) canonical release lives on
  # GitHub Landrylab/FungAMR (kept up to date). Current data tag: 070425 (2025-04-07).
  # NB: no bundled reference FASTA — ortholog alignments are in the paper's Suppl.
  # Data 1; the searchable copy is on CARD (card.mcmaster.ca/fungamrhome).
  # MARDy (mardy.net) has NO anonymous bulk export: its REST API (/api/) requires an
  # X-API-Key, and its content is largely subsumed by FungAMR — so it is not fetched here.
  local base="https://raw.githubusercontent.com/Landrylab/FungAMR/main"
  local ok=1
  ( cd "$DB/fungamr"
    for f in FungAMR_070425.tsv FungAMR_070425.csv FungAMR_070425.xlsx \
             drugs_class.csv genes_class.csv species_class.csv \
             curation_sheet_FungAMR.xlsx README.md; do
      dl "$base/${f// /%20}" "$f" || exit 1
    done
  ) >>"$LOGDIR/fungamr.log" 2>&1 || ok=0
  if [ "$ok" != 1 ]; then fail fungamr "FungAMR github download"; return; fi
  # FungAMR ships no sequences — build the alignment references af_resistance.py needs by
  # fetching each ref_seq_uniprot_accession from UniProt -> fungamr/reference_proteins.faa.
  log "building FungAMR reference proteins from UniProt (build_fungamr_refs.py)"
  if python3 "$(dirname "$0")/build_fungamr_refs.py" --data-dir "$DB" >>"$LOGDIR/fungamr.log" 2>&1 \
     && [ -s "$DB/fungamr/reference_proteins.faa" ]; then
    mark fungamr "tables + $(grep -c '>' "$DB/fungamr/reference_proteins.faa") reference proteins (MARDy=api-key-gated, skipped)"
  else
    fail fungamr "FungAMR reference-protein build (UniProt) failed — tables present, re-run to retry"
  fi
}

# ---- RVDB / mycoviral RdRp for mycovirus screen (Stage 10) -------------------
step_rvdb(){
  is_done rvdb && { log "rvdb present — skip"; return; }
  # Default = RVDB-prot v31.0 (2026-01) clustered/unique protein FASTA (100% id,
  # ~786k seqs, xz). Full unclustered set: U-RVDBv31.0-prot.fasta.xz. Check
  # https://rvdb-prot.pasteur.fr/ for newer versions.
  local RVDB_URL="${RVDB_URL:-https://rvdb-prot.pasteur.fr/files/U-RVDBv31.0-prot_unique.fasta.xz}"
  if [ -z "$RVDB_URL" ]; then fail rvdb "RVDB_URL not set — resolve current RVDB-prot release"; return; fi
  ( cd "$DB/rvdb" && dl "$RVDB_URL" rvdb.fasta.xz ) >>"$LOGDIR/rvdb.log" 2>&1 \
    && mark rvdb "$RVDB_URL" || fail rvdb "$RVDB_URL"
}

# ---- reference genomes for the assembly benchmarks (W1.1; small) -------------
# A. fumigatus A1163 (CEA10 lineage) for the CEA10 hybrid benchmark and A. flavus NRRL 3357
# for the Illumina-only benchmark. NCBI Datasets zip -> <acc>.fna; bin/benchmark_assembly.py.
step_benchmarks(){
  is_done benchmarks && { log "benchmark references present — skip"; return; }
  local ok=1 acc
  mkdir -p "$DB/benchmarks"
  for acc in GCA_000150145.1 GCA_009017415.1; do
    [ -s "$DB/benchmarks/$acc.fna" ] && { log "have $acc — skip"; continue; }
    log "downloading $acc -> $DB/benchmarks/$acc.fna"
    ( cd "$DB/benchmarks" \
      && dl "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA&filename=$acc.zip" "$acc.zip" \
      && unzip -o -q "$acc.zip" -d "$acc.tmp" \
      && cat "$acc.tmp"/ncbi_dataset/data/$acc/*.fna > "$acc.fna" && rm -rf "$acc.tmp" "$acc.zip" ) >>"$LOGDIR/benchmarks.log" 2>&1 \
      || { fail benchmarks "$acc"; ok=0; }
  done
  [ "$ok" = 1 ] && mark benchmarks "A1163 GCA_000150145.1 + NRRL3357 GCA_009017415.1"
}

# ---- type-material reference sets for the secondary ID loci (W2.3) ------------
# NCBI records flagged "sequence from type" for CaM, BenA, TEF1, RPB2 and LSU (D1/D2); the
# per-locus BLAST databases are built inside the identification task, so no BLAST is needed here.
step_markers(){
  is_done markers && { log "marker reference sets present — skip"; return; }
  log "fetching type-material marker sets -> $DB/markers"
  python3 "$REPO_DIR/bin/fetch_marker_refs.py" --out-dir "$DB/markers" >>"$LOGDIR/markers.log" 2>&1 \
    && mark markers "NCBI type-material CaM/BenA/TEF1/RPB2/LSU $(date +%F)" || fail markers "fetch_marker_refs.py (see logs/markers.log)"
}

# ---- fungal PubMLST schemes for mlst (W2.3) ------------------------------------
step_mlst(){
  is_done mlst && { log "mlst schemes present — skip"; return; }
  log "fetching PubMLST fungal schemes -> $DB/mlst/pubmlst"
  python3 "$REPO_DIR/bin/fetch_mlst_schemes.py" --out-dir "$DB/mlst" >>"$LOGDIR/mlst.log" 2>&1 \
    && mark mlst "PubMLST afumigatus calbicans cglabrata ctropicalis ckrusei $(date +%F)" || fail mlst "fetch_mlst_schemes.py (see logs/mlst.log)"
}

# ---- W2.6 extras: dbCAN HMMs, PHI-base, EffectorP 3 --------------------------------
# dbCAN CAZyme family HMMs (the HMMER module of run_dbcan; searched with hmmsearch in the base image).
step_dbcan(){
  is_done dbcan && { log "dbcan HMMs present — skip"; return; }
  # the dbCAN site serves its files through a download script (plain paths return an HTML page)
  local ver="${DBCAN_VERSION:-V14}"
  local url="${DBCAN_URL:-https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-$ver.txt}"
  log "downloading dbCAN HMM database $ver -> $DB/dbcan"
  ( cd "$DB/dbcan" && dl "$url" "dbCAN-HMMdb.txt" && grep -q '^HMMER3' dbCAN-HMMdb.txt && hmmpress -f dbCAN-HMMdb.txt >/dev/null 2>&1 || true ) >>"$LOGDIR/dbcan.log" 2>&1 \
    && [ -s "$DB/dbcan/dbCAN-HMMdb.txt" ] && grep -q '^HMMER3' "$DB/dbcan/dbCAN-HMMdb.txt" && mark dbcan "dbCAN-HMMdb-$ver ($url)" || fail dbcan "$url"
}

# PHI-base pathogen–host interaction proteins (CC BY; FASTA with phenotype in the header).
step_phibase(){
  is_done phibase && { log "phibase present — skip"; return; }
  local url="${PHIBASE_URL:-https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.fas}"
  log "downloading PHI-base -> $DB/phibase"
  ( cd "$DB/phibase" && dl "$url" "phi-base_current.fas" ) >>"$LOGDIR/phibase.log" 2>&1 \
    && [ -s "$DB/phibase/phi-base_current.fas" ] && mark phibase "$url" || fail phibase "$url"
}

# EffectorP 3.0 (GPL; Python + bundled WEKA, run with the Java of the base image).
step_effectorp(){
  is_done effectorp && { log "effectorp present — skip"; return; }
  local repo="${EFFECTORP_REPO:-https://github.com/JanaSperschneider/EffectorP-3.0}"
  log "cloning EffectorP 3.0 -> $DB/effectorp"
  ( rm -rf "$DB/effectorp/EffectorP-3.0" && git clone -q --depth 1 "$repo" "$DB/effectorp/EffectorP-3.0" \
    && cd "$DB/effectorp/EffectorP-3.0" && unzip -o -q weka-3-8-4.zip ) >>"$LOGDIR/effectorp.log" 2>&1 \
    && [ -f "$DB/effectorp/EffectorP-3.0/EffectorP.py" ] && [ -d "$DB/effectorp/EffectorP-3.0/weka-3-8-4" ] \
    && mark effectorp "$repo ($(git -C "$DB/effectorp/EffectorP-3.0" rev-parse --short HEAD))" || fail effectorp "$repo"
}

# ---- geNomad database (W2.8, stage 10) ---------------------------------------
# Zenodo release matching the pinned genomad (1.12 uses database v1.9); extracted to $DB/genomad_db.
step_genomad(){
  is_done genomad_db && { log "genomad_db present — skip"; return; }
  local url="${GENOMAD_DB_URL:-https://zenodo.org/api/records/14886553/files/genomad_db_v1.9.tar.gz/content}"
  mkdir -p "$DB/genomad_db"
  log "downloading geNomad database v1.9 (~0.8 GB) -> $DB/genomad_db"
  ( cd "$DB" && dl "$url" genomad_db.tar.gz && verify_archive genomad_db.tar.gz && tar xzf genomad_db.tar.gz && rm -f genomad_db.tar.gz ) >>"$LOGDIR/genomad.log" 2>&1 \
    && [ -s "$DB/genomad_db/genomad_db" ] || [ -s "$DB/genomad_db/version.txt" ] && mark genomad_db "$url" || fail genomad_db "$url"
}

# ---- reference genome set for genome-ANI novelty (W2.9, stage 12) -----------------
# The reference genome of every species in the genera of fungiforge/resources/novelty_genera.txt
# (NCBI Datasets; override with GENOMES_GENERA="Aspergillus,Candida" or add GENOMES_ACCESSIONS).
# Tens of GB for the default list; resumable (present files are kept).
step_genomes(){
  is_done refseq_fungi_genomes && { log "reference genome set present — skip"; return; }
  mkdir -p "$DB/refseq_fungi_genomes"
  log "fetching the reference genome set -> $DB/refseq_fungi_genomes"
  python3 "$REPO_DIR/bin/fetch_reference_genomes.py" --out-dir "$DB/refseq_fungi_genomes" \
      ${GENOMES_GENERA:+--genera "$GENOMES_GENERA"} --genera-file "$REPO_DIR/fungiforge/resources/novelty_genera.txt" \
      ${GENOMES_ACCESSIONS:+--accessions "$GENOMES_ACCESSIONS"} >>"$LOGDIR/genomes.log" 2>&1 \
    && mark refseq_fungi_genomes "$(ls "$DB/refseq_fungi_genomes"/*.fna 2>/dev/null | wc -l | tr -d ' ') genomes (fetch_reference_genomes.py)" \
    || fail refseq_fungi_genomes "fetch_reference_genomes.py (see logs/genomes.log)"
}

# ---- driver -----------------------------------------------------------------
STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=(images antismash funannotate eggnog busco unite kraken2 refseq_fungi fungamr rvdb benchmarks markers mlst dbcan phibase effectorp genomad)   # genomes (tens of GB) and interproscan: on request (6.9 GB + hours per genome)
log "==== fungiforge fetch_references start · DB=$DB · steps: ${STEPS[*]} ===="
for s in "${STEPS[@]}"; do "step_$s"; done
log "==== fetch_references finished ===="
# CURRENT state, one row per database: $MANIFEST is an append-only log, so a database that failed
# in July and succeeded in September has both rows in it and the raw log reads as broken. The
# `.done` marker on disk is the authority; the log stays for the audit trail.
log "current state (latest entry per database; full log: $MANIFEST)"
{
  echo -e "database\tstate\tsize\tdetail\twhen"
  for d in "$DB"/*/; do
    name="$(basename "$d")"
    [ "$name" = "logs" ] && continue
    last="$(awk -F'\t' -v n="$name" '$1==n {row=$0} END{print row}' "$MANIFEST")"
    state="missing"; [ -f "$d/.done" ] && state="present"
    [ -z "$(ls -A "$d" 2>/dev/null)" ] && state="empty"
    # these are staged only when asked for (large, and nothing needs them by default)
    case "$name" in
      interproscan|refseq_fungi_genomes)
        [ "$state" = "present" ] || state="on request" ;;
    esac
    size="$(du -sh "$d" 2>/dev/null | cut -f1)"
    echo -e "$name\t$state\t${size:--}\t$(echo "$last" | cut -f2 | cut -c1-60)\t$(echo "$last" | cut -f4)"
  done
} | column -t -s $'\t' 2>/dev/null || cat "$MANIFEST"
log "'missing'/'empty' -> re-fetch by naming the step: bin/fetch_references.sh <step>"
log "'on request' -> optional, nothing needs it by default: interproscan (6.9 GB, stage 07c), genomes (tens of GB, stage 12 ANI novelty)"
