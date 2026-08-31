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
#        fungamr rvdb
#
# NB several DBs are downloaded THROUGH their tool container (funannotate setup,
# eggnog, antismash) so the relevant image is pulled first. Version/URL-sensitive
# DBs (UNITE, Kraken2, sourmash-fungi) resolve the CURRENT release at run time;
# if a URL 404s the step logs and skips rather than guessing.
# ==============================================================================
set -uo pipefail

DB="${FUNGIFORGE_DB:?set FUNGIFORGE_DB to the target data dir, e.g. /Volumes/Extreme\ SSD/fungiforge_db}"
mkdir -p "$DB"/{containers,funannotate,eggnog,antismash,busco,unite,kraken2,refseq_fungi,fungamr,rvdb,logs}
LOGDIR="$DB/logs"; MANIFEST="$DB/MANIFEST.tsv"
[ -f "$MANIFEST" ] || echo -e "database\tdetail\tstatus\ttimestamp" > "$MANIFEST"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/fetch.log" ; }

# Robust, resumable, multi-connection download. aria2c preferred (16 parallel
# connections, --continue resumes a broken transfer from where it stopped, infinite
# retries); curl -C - is the fallback. Never restarts a partial file from scratch.
dl(){  # dl <url> <output-filename>   (run inside the target dir)
  local url="$1" out="$2"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c --continue=true -x16 -s16 -k1M --max-tries=0 --retry-wait=10 \
           --file-allocation=none --auto-file-renaming=false --console-log-level=warn \
           -o "$out" "$url"
  else
    curl -fL -C - --retry 999 --retry-delay 10 --retry-all-errors -o "$out" "$url"
  fi
}

is_done(){ [ -f "$DB/$1/.done" ]; }
mark(){ date '+%F %T' > "$DB/$1/.done"; echo -e "$1\t$2\tOK\t$(date '+%F %T')" >> "$MANIFEST"; log "✓ $1 done ($2)"; }
fail(){ echo -e "$1\t$2\tFAILED\t$(date '+%F %T')" >> "$MANIFEST"; log "✗ $1 FAILED ($2) — continuing"; }

# ---- container images (linux/amd64; emulated on Mac, native on HPC) ----------
IMAGES=(
  "staphb/flye:latest"
  "staphb/medaka:latest"
  "staphb/filtlong:latest"
  "antismash/standalone:8.0.0"
  "nextgenusfs/funannotate:latest"      # Funannotate: predict + annotate + setup
  "dfam/tetools:latest"                 # RepeatModeler2 + RepeatMasker
  "ezlabgva/busco:v5.7.1_cv1"           # BUSCO (fallback if compleasm unavailable)
)
step_images(){
  is_done containers && { log "images already pulled — skip"; return; }
  local ok=1
  for img in "${IMAGES[@]}"; do
    log "docker pull $img"
    docker pull --platform linux/amd64 "$img" >>"$LOGDIR/images.log" 2>&1 || { fail containers "$img"; ok=0; }
  done
  [ "$ok" = 1 ] && mark containers "${#IMAGES[@]} images"
}

# ---- antiSMASH databases (~9 GB) via the antismash 8 image -------------------
step_antismash(){
  is_done antismash && { log "antismash db present — skip"; return; }
  log "downloading antiSMASH databases -> $DB/antismash"
  # The image entrypoint is `antismash`; the DB downloader is a separate console
  # script, so override the entrypoint. `download-antismash-databases` takes --database-dir.
  docker run --rm --platform linux/amd64 --entrypoint download-antismash-databases \
      -v "$DB/antismash":/db antismash/standalone:8.0.0 --database-dir /db >>"$LOGDIR/antismash.log" 2>&1 \
    && mark antismash "antismash8 db" || fail antismash "download-antismash-databases"
}

# ---- Funannotate database (~30-50 GB): Pfam, dbCAN, MEROPS, InterPro, BUSCO --
step_funannotate(){
  is_done funannotate && { log "funannotate db present — skip"; return; }
  log "funannotate setup -i all -> $DB/funannotate (large, hours)"
  docker run --rm --platform linux/amd64 -e FUNANNOTATE_DB=/data -v "$DB/funannotate":/data \
      nextgenusfs/funannotate:latest funannotate setup -i all -d /data >>"$LOGDIR/funannotate.log" 2>&1 \
    && mark funannotate "funannotate setup all" || fail funannotate "funannotate setup"
}

# ---- eggNOG database (~50 GB) -----------------------------------------------
# eggnog-mapper ships download_eggnog_data.py; run it from the funannotate image
# (bundles eggnog-mapper) or the eggnog-mapper image. -y = accept all downloads.
step_eggnog(){
  is_done eggnog && { log "eggnog db present — skip"; return; }
  log "downloading eggNOG data -> $DB/eggnog (large)"
  docker run --rm --platform linux/amd64 -v "$DB/eggnog":/eggnog \
      nextgenusfs/funannotate:latest download_eggnog_data.py -y --data_dir /eggnog >>"$LOGDIR/eggnog.log" 2>&1 \
    && mark eggnog "eggnog data" || fail eggnog "download_eggnog_data.py"
}

# ---- BUSCO / compleasm fungal lineages --------------------------------------
# compleasm is fastest; if not installed, fall back to BUSCO container download.
step_busco(){
  is_done busco && { log "busco lineages present — skip"; return; }
  log "fetching fungal BUSCO/compleasm lineages -> $DB/busco"
  if command -v compleasm >/dev/null 2>&1; then
    compleasm download fungi_odb10 -L "$DB/busco" >>"$LOGDIR/busco.log" 2>&1 \
      && mark busco "compleasm fungi_odb10" || fail busco "compleasm download"
  else
    docker run --rm --platform linux/amd64 -v "$DB/busco":/busco -w /busco \
        ezlabgva/busco:v5.7.1_cv1 busco --download fungi_odb10 >>"$LOGDIR/busco.log" 2>&1 \
      && mark busco "busco fungi_odb10" || fail busco "busco --download"
  fi
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
  ( cd "$DB/unite" && dl "$UNITE_URL" unite.tgz && tar xzf unite.tgz ) >>"$LOGDIR/unite.log" 2>&1 \
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
  ( cd "$DB/kraken2" && dl "$K2_URL" k2.tgz && tar xzf k2.tgz ) >>"$LOGDIR/kraken2.log" 2>&1 \
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

# ---- driver -----------------------------------------------------------------
STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=(images antismash funannotate eggnog busco unite kraken2 refseq_fungi fungamr rvdb)
log "==== fungiforge fetch_references start · DB=$DB · steps: ${STEPS[*]} ===="
for s in "${STEPS[@]}"; do "step_$s"; done
log "==== fetch_references finished ===="
log "manifest: $MANIFEST"; column -t -s $'\t' "$MANIFEST" 2>/dev/null || cat "$MANIFEST"
