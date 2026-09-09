#!/usr/bin/env bash
###############################################################################
# fungiforge preflight_qc.sh
#
# WHAT THIS IS FOR
#   A fast, read-only triage you run on FRESH Oxford Nanopore (ONT) data BEFORE
#   committing to a long fungiforge run. On this Apple-Silicon Mac the pipeline
#   is emulated (linux/amd64) and a real run is multi-hour, so you want to know
#   up front, per sample:
#     - Is there enough usable data to assemble a genome? (GO / MARGINAL / NO-GO)
#     - Is this actually whole-genome shotgun, or an ITS amplicon run that
#       fungiforge is the WRONG tool for?
#     - WHAT is each isolate? (a read-based species ID that works even when
#       coverage is too low for assembly to succeed — the "salvage" step.)
#
# QUICK START (brand-new dataset, per-barcode ONT chunks under one dir):
#   export FUNGIFORGE_DB=/path/to/fungiforge_db
#   bin/preflight_qc.sh \
#       --input "/path/to/Passed samples" \
#       --merged-out ./preflight_merged \
#       --outdir ./preflight_qc \
#       --genome-size 30000000
#   # (input path may contain spaces; --merged-out / --outdir must NOT)
#
# THE 8 CHECKS (each section below prints its own interpretation):
#   1. ENVIRONMENT READINESS  — docker image present; DB root + required subdirs
#      (unite/ funannotate/ antismash/ busco/lineages/ fungamr/reference_proteins.faa);
#      free disk; LOUD warning if the DB or work/read path contains a SPACE
#      (Flye/funannotate/antiSMASH break on spaces); count of macOS AppleDouble
#      "._*" files in the DB (they break antiSMASH/funannotate).
#   2. MERGE PER-BARCODE CHUNKS — ONT arrives as many *.fastq.gz per barcode dir.
#      Concatenate each barcode dir into one <SAMPLE>.ont.fastq.gz in a SPACE-FREE
#      output dir (gzip members concatenate with plain cat). Default name map:
#      barcode74 -> BC74 (see --sample-prefix). Use --skip-merge if already merged.
#   3. READ STATS — seqkit stats -a: num_seqs, sum_len, avg_len, read N50, Q20%,
#      Q30%, GC%. The headline "how much / how good" table.
#   4. READ-LENGTH BUCKETS — counts in <300 / 300-500 / 500-800 / 800-1200 /
#      1200-3000 / >3000 bp. Eyeball the SHAPE: a spike in one narrow band =>
#      amplicon-like; a broad smear to long reads => shotgun; heavy <300 bp =>
#      short-read burden / junk that won't assemble.
#   5. AMPLICON-vs-WGS DISCRIMINATOR — subsample ~2000 reads (>=200bp), vsearch
#      against UNITE at --id 0.75. Report the FRACTION that hit ITS.
#        * HIGH hit fraction (>50%), concentrated on ONE species => ITS AMPLICON
#          run: fungiforge is the WRONG tool — use an ITS/metabarcoding workflow.
#        * VERY LOW fraction (<~5%), scattered across many taxa => whole-genome
#          shotgun => correct input for fungiforge.
#   6. ASSEMBLY-USABLE YIELD & ROUGH COVERAGE — reads/bp at >=1000 and >=3000 bp,
#      and fold-coverage for --genome-size (default 30 Mb, molds; yeasts ~15 Mb so
#      DOUBLE the number). Verdict from usable (>=1kb) coverage. HEURISTIC guide:
#          >=30x usable  => GO
#          10-30x usable => MARGINAL (expect a fragmented assembly)
#          <10x usable   => NO-GO (assembly likely fails; rely on check 7)
#   7. READ-BASED SPECIES ID (assembly-free, the SALVAGE step) — reads >=300bp,
#      vsearch against UNITE at --id 0.85, tally the top target species. Prints the
#      top 3 species + hit counts per sample. Works even at low coverage because
#      fungal rDNA/ITS is tandemly repeated ~50-200x, so ITS-bearing reads are
#      enriched and present long before genome coverage is assembly-worthy.
#   8. ARTIFACTS — writes a per-sample summary TSV and a ready-to-edit fungiforge
#      samplesheet stub (ont_fastq -> merged files; metadata columns = NA).
#
# TOOLING FACTS (tested on this machine — do not "simplify" these away):
#   * Everything runs inside the base image aduton1000/fungiforge:0.1.0.
#   * arm64 host + linux/amd64 image => every docker run passes
#     --platform linux/amd64 and -e HOME=/tmp -e USER=fungiforge.
#   * Binaries live at /opt/conda/bin/<tool> (NOT on the bash -lc PATH) -> full paths.
#   * kraken2 is NOT in the image -> read-based ID uses vsearch vs UNITE.
#   * bc is NOT in the image -> all arithmetic is done in awk.
#   * Globs like *.fastq.gz are expanded INSIDE the container.
#
#   Read this header with:  head -80 bin/preflight_qc.sh
#   Full walkthrough:        docs/preflight_qc.md
###############################################################################
set -euo pipefail

# ----------------------------------------------------------------------------
# Defaults
# ----------------------------------------------------------------------------
IMAGE="${FUNGIFORGE_IMAGE:-aduton1000/fungiforge:0.1.0}"
SIF="${FUNGIFORGE_SIF:-}"   # set (or --sif) => run via apptainer/singularity instead of docker (HPC)
DB="${FUNGIFORGE_DB:-}"
INPUT=""
MERGED="./preflight_merged"
OUTDIR="./preflight_qc"
GENOME_SIZE=30000000
SAMPLE_PREFIX="BC"
SKIP_MERGE=0
FORCE=0
AMPLICON_SUBN=2000     # check 5: reads subsampled for the amplicon-vs-WGS test
ID_MAX_READS=50000     # check 7: cap reads fed to vsearch to keep ID fast

# ---- report/verdict thresholds (kept explicit so the report can cite them) ----
COV_GO=30              # >=1kb fold-cov at/above this => GO for de novo assembly
COV_MARGINAL=10        # >=1kb fold-cov in [COV_MARGINAL,COV_GO) => MARGINAL
GOOD_N50=5000          # read N50 at/above this = "good ONT"; below => flag
AMP_PCT_HI=50          # check 5: %ITS above this + few species => amplicon
AMP_PCT_LO=5           # check 5: %ITS below this => shotgun

# ----------------------------------------------------------------------------
# Usage
# ----------------------------------------------------------------------------
usage() {
  cat <<'EOF'
preflight_qc.sh — fast ONT triage before a long fungiforge run

USAGE
  bin/preflight_qc.sh --input DIR [options]

REQUIRED
  --input DIR          Input reads. Either (a) a dir containing per-barcode
                       subdirs each holding *.fastq.gz chunks (default; will be
                       merged), or (b) a dir of already-merged per-sample
                       *.fastq.gz with --skip-merge. The input path MAY contain
                       spaces (e.g. ".../Passed samples").
  --db DIR             fungiforge DB root (or set env FUNGIFORGE_DB). Must
                       contain unite/, funannotate/, antismash/, busco/lineages/,
                       fungamr/reference_proteins.faa.

OPTIONS
  --merged-out DIR     Where merged <SAMPLE>.ont.fastq.gz are written.
                       MUST be space-free. Default: ./preflight_merged
  --outdir DIR         Where summary TSV + samplesheet stub go.
                       Default: ./preflight_qc
  --genome-size N      Expected genome size (bp) for coverage. Default 30000000
                       (molds). Yeasts ~15000000 — halve this or double the cov.
  --sample-prefix STR  Prefix for barcodeNN -> <STR>NN naming. Default: BC
  --skip-merge         Input already holds merged per-sample *.fastq.gz; skip
                       the merge step and read them directly.
  --amplicon-subn N    Reads subsampled for the amplicon test. Default 2000.
  --id-max-reads N     Cap on reads fed to the species-ID vsearch. Default 50000.
  --image NAME         Base docker image. Default aduton1000/fungiforge:0.1.0
  --sif PATH           Run inside this .sif with apptainer/singularity instead of
                       docker (HPC). Default $FUNGIFORGE_SIF.
  --force              Re-merge / overwrite existing merged outputs.
  -h, --help           This help.

QUICK START
  export FUNGIFORGE_DB=/path/to/fungiforge_db
  bin/preflight_qc.sh --input "/path/to/Passed samples" \
      --merged-out ./preflight_merged --outdir ./preflight_qc

Full documentation: docs/preflight_qc.md
EOF
}

# ----------------------------------------------------------------------------
# Arg parsing
# ----------------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --input)         INPUT="${2:-}"; shift 2 ;;
    --db)            DB="${2:-}"; shift 2 ;;
    --merged-out)    MERGED="${2:-}"; shift 2 ;;
    --outdir)        OUTDIR="${2:-}"; shift 2 ;;
    --genome-size)   GENOME_SIZE="${2:-}"; shift 2 ;;
    --sample-prefix) SAMPLE_PREFIX="${2:-}"; shift 2 ;;
    --amplicon-subn) AMPLICON_SUBN="${2:-}"; shift 2 ;;
    --id-max-reads)  ID_MAX_READS="${2:-}"; shift 2 ;;
    --image)         IMAGE="${2:-}"; shift 2 ;;
    --sif)           SIF="${2:-}"; shift 2 ;;
    --skip-merge)    SKIP_MERGE=1; shift ;;
    --force)         FORCE=1; shift ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; echo "Run with --help." >&2; exit 2 ;;
  esac
done

# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
hr()   { printf '%s\n' "======================================================================"; }
banner(){ hr; echo "$1"; hr; }

has_space(){ case "$1" in *' '*) return 0 ;; *) return 1 ;; esac; }

abspath(){
  # portable absolute path (macOS lacks GNU realpath by default)
  local p="$1"
  if [ -d "$p" ]; then ( cd "$p" && pwd )
  else ( cd "$(dirname "$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")" ); fi
}

# map a barcode dir name to a sample name: barcode74 -> BC74; else sanitized name
map_name(){
  local b="$1"
  if [[ "$b" =~ ^barcode([0-9]+)$ ]]; then
    printf '%s%s\n' "$SAMPLE_PREFIX" "${BASH_REMATCH[1]}"
  else
    printf '%s\n' "${b// /_}"
  fi
}

# pull a named column value out of a 2-line (header+data) seqkit TSV block
col(){ # $1 = block text, $2 = column name
  printf '%s\n' "$1" | awk -F'\t' -v c="$2" '
    NR==1{for(i=1;i<=NF;i++) if($i==c) k=i}
    NR==2{print (k?$k:"NA")}'
}
num(){ local v="${1:-}"; [ -z "$v" ] && echo 0 || echo "$v"; }

# ----------------------------------------------------------------------------
# docker wrapper — one place for the platform/env/mount incantation.
# Forwards BASE/UNITE/AMPLICON_SUBN/ID_MAX_READS from the host env so inner
# scripts stay single-quoted (no fragile nested quoting).
# ----------------------------------------------------------------------------
dock(){ # $1 = inner bash script (single-quoted at call site)
  if [ -n "$SIF" ]; then
    # apptainer/singularity: same /reads,/db,/out layout; host env is inherited so the
    # -e forwards are implicit. --containall keeps $HOME/CWD out, mirroring docker.
    export BASE UNITE AMPLICON_SUBN ID_MAX_READS COV_GO GOOD_N50 AMP_PCT_HI
    "$CRT" exec --containall --env USER=fungiforge --pwd /out \
      -B "$MERGED_ABS":/reads:ro -B "$DB":/db:ro -B "$OUTDIR_ABS":/out \
      "$SIF" bash -c "$1"
  else
    docker run --rm --platform linux/amd64 \
      -e HOME=/tmp -e USER=fungiforge \
      -e BASE -e UNITE -e AMPLICON_SUBN -e ID_MAX_READS \
      -e COV_GO -e GOOD_N50 -e AMP_PCT_HI \
      -v "$MERGED_ABS":/reads:ro \
      -v "$DB":/db:ro \
      -v "$OUTDIR_ABS":/out \
      "$IMAGE" bash -c "$1"
  fi
}

# ----------------------------------------------------------------------------
# Report accumulators (filled as checks run; consumed by the report section)
# ----------------------------------------------------------------------------
declare -a ENV_ROWS=()          # "STATUS|label|detail"  STATUS in PASS/WARN/FAIL
envrow(){ ENV_ROWS+=("$1|$2|$3"); }
declare -a R_SAMPLE R_NALL R_BPALL R_N50 R_GC R_PITS R_COVALL R_COV1K \
           R_VERDICT R_AMP R_TOP R_REASON R_B=()
RUN_DATE="$(date '+%Y-%m-%d %H:%M:%S %Z')"
RUN_HOST="$(hostname 2>/dev/null || echo unknown)"

###############################################################################
# CHECK 1 — ENVIRONMENT READINESS
###############################################################################
banner "CHECK 1  ENVIRONMENT READINESS"

CRT=""
if [ -n "$SIF" ]; then
  command -v apptainer >/dev/null 2>&1 && CRT=apptainer
  [ -z "$CRT" ] && command -v singularity >/dev/null 2>&1 && CRT=singularity
  [ -n "$CRT" ] || die "--sif given but neither apptainer nor singularity is on PATH."
  [ -s "$SIF" ] || die "container image not found: $SIF"
  echo "[ok] container image present ($CRT): $SIF"
  envrow PASS "container image" "$SIF ($CRT)"
else
  command -v docker >/dev/null 2>&1 || die "docker not found on PATH. Install/start Docker Desktop (or pass --sif on an HPC)."
  docker info >/dev/null 2>&1 || die "docker daemon not responding. Is Docker Desktop running?"

  # NB: capture first — piping `docker images` into `grep -q` makes grep close the
  # pipe on first match, which SIGPIPEs `docker images` and trips `pipefail`.
  IMG_REPO="${IMAGE%:*}"
  IMG_LIST="$(docker images --format '{{.Repository}}:{{.Tag}}' 2>/dev/null || true)"
  if printf '%s\n' "$IMG_LIST" | grep -q "$IMG_REPO"; then
    echo "[ok] base image present: $IMAGE"
    envrow PASS "docker image" "$IMAGE"
  else
    die "base image '$IMAGE' not found (docker images | grep fungiforge). Build/pull it first."
  fi
fi

[ -n "$INPUT" ] || die "--input is required. Run with --help."
[ -e "$INPUT" ] || die "--input path does not exist: $INPUT"
[ -n "$DB" ]    || die "--db not set and FUNGIFORGE_DB is empty. Point to the fungiforge DB root."
[ -d "$DB" ]    || die "--db path is not a directory: $DB"

echo "[..] checking DB layout under: $DB"
missing=0
for need in unite funannotate antismash busco/lineages fungamr/reference_proteins.faa; do
  if [ -e "$DB/$need" ]; then
    echo "     [ok] $need"
  else
    echo "     [MISSING] $need"; missing=1
  fi
done
[ "$missing" -eq 0 ] || die "DB is incomplete (see [MISSING] above). Run bin/fetch_references.sh."
envrow PASS "DB layout" "unite, funannotate, antismash, busco/lineages, fungamr present"

# Resolve the UNITE fasta (filename changes with DB version). Prefer the
# non-_dev general_release_dynamic file; fall back gracefully.
UNITE=""
for f in "$DB"/unite/sh_general_release_dynamic_*.fasta; do
  [ -e "$f" ] || continue
  case "$f" in *_dev.fasta) continue ;; esac
  UNITE="$(basename "$f")"; break
done
if [ -z "$UNITE" ]; then
  for f in "$DB"/unite/*.fasta; do [ -e "$f" ] && { UNITE="$(basename "$f")"; break; }; done
fi
[ -n "$UNITE" ] || die "no UNITE fasta found under $DB/unite/ (expected sh_general_release_dynamic_*.fasta)."
echo "[ok] UNITE reference: $UNITE"
envrow PASS "UNITE reference" "$UNITE"

# SPACE warnings — the documented Flye/funannotate/antiSMASH gotcha.
if has_space "$DB"; then
  warn "DB path contains a SPACE: '$DB'"
  warn "  -> Flye/funannotate/antiSMASH break on spaces. Move the DB to a space-free path."
  envrow WARN "DB path has space" "$DB"
else
  envrow PASS "DB path space-free" "$DB"
fi
if has_space "$MERGED"; then
  warn "--merged-out contains a SPACE: '$MERGED'"
  warn "  -> merged reads MUST live on a space-free path. Choose another --merged-out."
  envrow WARN "merged path has space" "$MERGED"
else
  envrow PASS "merged path space-free" "$MERGED"
fi
if [ "$SKIP_MERGE" -eq 1 ] && has_space "$INPUT"; then
  warn "--skip-merge with a spaced --input ('$INPUT'): the merged reads fungiforge"
  warn "  will consume already contain a space and will break the real run."
fi

# macOS AppleDouble files in the DB (break antiSMASH/funannotate)
adcount=$(find "$DB" -name '._*' 2>/dev/null | wc -l | tr -d ' ')
if [ "$adcount" -gt 0 ]; then
  warn "found $adcount macOS AppleDouble '._*' file(s) under the DB."
  warn "  -> these break antiSMASH/funannotate. Clean with:  dot_clean -m '$DB'  (or: find '$DB' -name '._*' -delete)"
  envrow WARN "AppleDouble '._*' in DB" "$adcount file(s) — dot_clean -m '$DB'"
else
  echo "[ok] no macOS '._*' AppleDouble files in the DB"
  envrow PASS "no AppleDouble '._*' in DB" "0 files"
fi

echo "[..] free disk:"
df -h "$DB" "${MERGED%/*}" 2>/dev/null | sed 's/^/     /' || df -h | sed 's/^/     /'
DISK_AVAIL="$(df -h "$DB" 2>/dev/null | awk 'NR==2{print $4}')"
envrow PASS "free disk (DB volume)" "${DISK_AVAIL:-unknown} available"

###############################################################################
# CHECK 2 — MERGE PER-BARCODE CHUNKS
###############################################################################
banner "CHECK 2  MERGE PER-BARCODE CHUNKS"

mkdir -p "$MERGED" "$OUTDIR"
MERGED_ABS="$(abspath "$MERGED")"
OUTDIR_ABS="$(abspath "$OUTDIR")"

if [ "$SKIP_MERGE" -eq 1 ]; then
  echo "[..] --skip-merge: reading already-merged fastqs directly from --input"
  # Point the merged dir at the input for the rest of the run.
  MERGED_ABS="$(abspath "$INPUT")"
  shopt -s nullglob
  merged_files=("$MERGED_ABS"/*.fastq.gz)
  shopt -u nullglob
  [ "${#merged_files[@]}" -gt 0 ] || die "--skip-merge set but no *.fastq.gz found in $INPUT"
  echo "[ok] found ${#merged_files[@]} merged fastq(s) in $INPUT"
else
  shopt -s nullglob
  had_any=0
  for d in "$INPUT"/*/; do
    [ -d "$d" ] || continue
    files=("$d"*.fastq.gz)
    [ "${#files[@]}" -gt 0 ] || continue
    had_any=1
    bn="$(basename "$d")"
    sample="$(map_name "$bn")"
    out="$MERGED_ABS/${sample}.ont.fastq.gz"
    if [ -e "$out" ] && [ "$FORCE" -eq 0 ]; then
      echo "     [skip] $bn -> $sample  (exists; use --force to redo)"
      continue
    fi
    echo "     [merge] $bn -> ${sample}.ont.fastq.gz  (${#files[@]} chunk(s))"
    cat "${files[@]}" > "$out"     # gzip members concatenate cleanly
  done
  shopt -u nullglob
  [ "$had_any" -eq 1 ] || die "no per-barcode subdirs with *.fastq.gz found under: $INPUT
  (If these are already merged per-sample fastqs, re-run with --skip-merge.)"
  echo "[ok] merged reads in: $MERGED_ABS"
fi

# Discover the merged samples we will report on.
shopt -s nullglob
SAMPLE_FILES=("$MERGED_ABS"/*.fastq.gz)
shopt -u nullglob
[ "${#SAMPLE_FILES[@]}" -gt 0 ] || die "no merged *.fastq.gz to analyse in $MERGED_ABS"
echo "[ok] ${#SAMPLE_FILES[@]} sample(s) to profile"

sample_of(){ # basename -> sample name
  local b; b="$(basename "$1")"
  local s="${b%.ont.fastq.gz}"; [ "$s" = "$b" ] && s="${b%.fastq.gz}"
  printf '%s\n' "$s"
}

# Stage awk helpers into the mounted outdir (kept out of quoting hell).
mkdir -p "$OUTDIR_ABS/.scripts"
cat > "$OUTDIR_ABS/.scripts/buckets.awk" <<'AWK'
# input: seqkit fx2tab -nl  ->  name<TAB>length
BEGIN{OFS="\t"}
{ l=$2
  if(l<300)        b1++
  else if(l<500)   b2++
  else if(l<800)   b3++
  else if(l<1200)  b4++
  else if(l<3000)  b5++
  else             b6++ }
END{ print b1+0,b2+0,b3+0,b4+0,b5+0,b6+0 }
AWK

# TSV accumulators
SUMMARY="$OUTDIR_ABS/preflight_summary.tsv"
printf 'sample\tnum_reads\tsum_bp\tread_N50\tGC\tpct_ITS_hits\tcov_all\tcov_ge1kb\tassembly_verdict\ttop_species\n' > "$SUMMARY"

###############################################################################
# Per-sample profiling loop (checks 3-7)
###############################################################################
for f in "${SAMPLE_FILES[@]}"; do
  BASE="$(basename "$f")"
  SAMPLE="$(sample_of "$f")"
  export BASE AMPLICON_SUBN ID_MAX_READS UNITE

  banner "SAMPLE  $SAMPLE   ($BASE)"

  # ---- CHECK 3 + 6 data: combined stats (ALL / >=1kb / >=3kb) ----
  STATS_BLOB="$(dock '
    set -e
    echo "@@ALL";  /opt/conda/bin/seqkit stats -a -T "/reads/$BASE"
    echo "@@GE1K"; /opt/conda/bin/seqkit seq -m 1000 "/reads/$BASE" 2>/dev/null | /opt/conda/bin/seqkit stats -T -
    echo "@@GE3K"; /opt/conda/bin/seqkit seq -m 3000 "/reads/$BASE" 2>/dev/null | /opt/conda/bin/seqkit stats -T -
  ')"

  get_block(){ printf '%s\n' "$STATS_BLOB" | awk -v m="@@$1" 'p&&c<2{print;c++} $0==m{p=1;c=0}'; }
  ALL_B="$(get_block ALL)"; GE1K_B="$(get_block GE1K)"; GE3K_B="$(get_block GE3K)"

  n_all=$(num "$(col "$ALL_B" num_seqs)")
  bp_all=$(num "$(col "$ALL_B" sum_len)")
  avg_all=$(num "$(col "$ALL_B" avg_len)")
  n50=$(num "$(col "$ALL_B" N50)")
  q20=$(num "$(col "$ALL_B" 'Q20(%)')")
  q30=$(num "$(col "$ALL_B" 'Q30(%)')")
  gc=$(num "$(col "$ALL_B" 'GC(%)')")
  n_1k=$(num "$(col "$GE1K_B" num_seqs)");  bp_1k=$(num "$(col "$GE1K_B" sum_len)")
  n_3k=$(num "$(col "$GE3K_B" num_seqs)");  bp_3k=$(num "$(col "$GE3K_B" sum_len)")

  echo
  echo "CHECK 3  READ STATS"
  { printf 'metric\tvalue\n'
    printf 'num_seqs\t%s\n' "$n_all"
    printf 'sum_len(bp)\t%s\n' "$bp_all"
    printf 'avg_len\t%s\n' "$avg_all"
    printf 'read_N50\t%s\n' "$n50"
    printf 'Q20(%%)\t%s\n' "$q20"
    printf 'Q30(%%)\t%s\n' "$q30"
    printf 'GC(%%)\t%s\n' "$gc"
  } | column -t -s "$(printf '\t')"

  # ---- CHECK 4: length-distribution buckets ----
  echo
  echo "CHECK 4  READ-LENGTH BUCKETS (bp)"
  BUCKETS="$(dock '/opt/conda/bin/seqkit fx2tab -nl "/reads/$BASE" 2>/dev/null | awk -f /out/.scripts/buckets.awk')"
  { printf '<300\t300-500\t500-800\t800-1200\t1200-3000\t>3000\n'
    printf '%s\n' "$BUCKETS"
  } | column -t -s "$(printf '\t')"
  echo "  (shape hint: one narrow spike => amplicon-like; broad smear to long reads => shotgun;"
  echo "   heavy <300 bp => short/junk that will not assemble)"

  # ---- CHECK 5: amplicon-vs-WGS discriminator ----
  echo
  echo "CHECK 5  AMPLICON-vs-WGS DISCRIMINATOR (subsample ~$AMPLICON_SUBN reads, UNITE id 0.75)"
  AMP="$(dock '
    /opt/conda/bin/seqkit seq -m 200 "/reads/$BASE" 2>/dev/null \
      | /opt/conda/bin/seqkit sample -n "$AMPLICON_SUBN" -s 11 2>/dev/null \
      | /opt/conda/bin/seqkit fq2fa 2>/dev/null > /tmp/sub.fa
    N=$(grep -c "^>" /tmp/sub.fa 2>/dev/null || echo 0)
    /opt/conda/bin/vsearch --usearch_global /tmp/sub.fa --db "/db/unite/$UNITE" \
      --id 0.75 --top_hits_only --strand both --maxaccepts 1 \
      --userout /tmp/amp.tsv --userfields query+target 2>/dev/null || true
    if [ -s /tmp/amp.tsv ]; then
      H=$(cut -f1 /tmp/amp.tsv | sort -u | wc -l | tr -d " ")
      NSP=$(cut -f2 /tmp/amp.tsv | cut -d"|" -f1 | sort -u | wc -l | tr -d " ")
      TOP=$(cut -f2 /tmp/amp.tsv | cut -d"|" -f1 | sort | uniq -c | sort -rnk1 | head -1 | awk "{print \$2\"(\"\$1\")\"}")
    else H=0; NSP=0; TOP="none"; fi
    echo "$N $H $NSP $TOP"
  ')"
  read -r amp_n amp_h amp_nsp amp_top <<<"$AMP"
  amp_n=$(num "$amp_n"); amp_h=$(num "$amp_h"); amp_nsp=$(num "$amp_nsp")
  pct_its=$(awk -v h="$amp_h" -v n="$amp_n" 'BEGIN{ printf "%.1f", (n>0? 100*h/n : 0) }')
  echo "  subsampled=$amp_n  ITS_hits=$amp_h  pct_ITS=${pct_its}%  distinct_species=$amp_nsp  top=$amp_top"
  amp_verdict=$(awk -v p="$pct_its" -v s="$amp_nsp" -v hi="$AMP_PCT_HI" -v lo="$AMP_PCT_LO" 'BEGIN{
      if(p>hi && s<=3) print "LIKELY ITS AMPLICON  -> fungiforge is the WRONG tool; use an ITS/metabarcoding workflow";
      else if(p<lo)    print "LIKELY WGS SHOTGUN   -> correct input for fungiforge";
      else             print "AMBIGUOUS            -> inspect buckets (check 4) + species spread (check 7)"}')
  amp_class=$(awk -v p="$pct_its" -v s="$amp_nsp" -v hi="$AMP_PCT_HI" -v lo="$AMP_PCT_LO" 'BEGIN{
      if(p>hi && s<=3) print "AMPLICON"; else if(p<lo) print "WGS"; else print "AMBIGUOUS"}')
  echo "  => $amp_verdict"

  # ---- CHECK 6: assembly-usable yield & rough coverage ----
  echo
  echo "CHECK 6  ASSEMBLY-USABLE YIELD & ROUGH COVERAGE (genome=${GENOME_SIZE} bp)"
  cov_all=$(awk -v s="$bp_all" -v g="$GENOME_SIZE" 'BEGIN{printf "%.1f",(g>0?s/g:0)}')
  cov_1k=$(awk  -v s="$bp_1k"  -v g="$GENOME_SIZE" 'BEGIN{printf "%.1f",(g>0?s/g:0)}')
  cov_3k=$(awk  -v s="$bp_3k"  -v g="$GENOME_SIZE" 'BEGIN{printf "%.1f",(g>0?s/g:0)}')
  { printf 'threshold\treads\tbp\tfold_cov\n'
    printf 'all\t%s\t%s\t%sx\n'    "$n_all" "$bp_all" "$cov_all"
    printf '>=1kb\t%s\t%s\t%sx\n'  "$n_1k"  "$bp_1k"  "$cov_1k"
    printf '>=3kb\t%s\t%s\t%sx\n'  "$n_3k"  "$bp_3k"  "$cov_3k"
  } | column -t -s "$(printf '\t')"
  verdict=$(awk -v c="$cov_1k" -v go="$COV_GO" -v mg="$COV_MARGINAL" 'BEGIN{ if(c>=go)print"GO"; else if(c>=mg)print"MARGINAL"; else print"NO-GO"}')
  echo "  verdict (from >=1kb coverage, HEURISTIC): $verdict"
  echo "  guide: >=30x GO | 10-30x MARGINAL (fragmented) | <10x NO-GO.  Yeasts ~15Mb => double the x."

  # ---- CHECK 7: read-based species ID (salvage) ----
  echo
  echo "CHECK 7  READ-BASED SPECIES ID (reads >=300bp, UNITE id 0.85; top 3)"
  IDOUT="$(dock '
    /opt/conda/bin/seqkit seq -m 300 "/reads/$BASE" 2>/dev/null \
      | /opt/conda/bin/seqkit sample -n "$ID_MAX_READS" -s 13 2>/dev/null \
      | /opt/conda/bin/seqkit fq2fa 2>/dev/null > /tmp/q.fa
    /opt/conda/bin/vsearch --usearch_global /tmp/q.fa --db "/db/unite/$UNITE" \
      --id 0.85 --top_hits_only --maxaccepts 3 --strand both \
      --userout /tmp/id.tsv --userfields target 2>/dev/null || true
    if [ -s /tmp/id.tsv ]; then
      cut -f1 /tmp/id.tsv | cut -d"|" -f1 | sort | uniq -c | sort -rnk1 | head -3
    else
      echo "      0 no_hit"
    fi
  ')"
  printf '%s\n' "$IDOUT" | awk '{cnt=$1; $1=""; sub(/^ +/,""); printf "  %6d  %s\n", cnt, $0}'
  top_species=$(printf '%s\n' "$IDOUT" | head -1 | awk '{cnt=$1;$1="";sub(/^ +/,"");printf "%s(%d)",$0,cnt}')
  [ -n "$top_species" ] || top_species="no_hit"

  # ---- per-sample verdict reason (thresholds explicit) ----
  reason=$(awk -v v="$verdict" -v c="$cov_1k" -v n="$n50" -v go="$COV_GO" -v mg="$COV_MARGINAL" -v gn="$GOOD_N50" 'BEGIN{
      if(v=="GO")       r=sprintf(">=1kb coverage %sx >= %dx", c, go);
      else if(v=="MARGINAL") r=sprintf(">=1kb coverage %sx in [%d,%d)x -> expect a fragmented assembly", c, mg, go);
      else              r=sprintf(">=1kb coverage %sx < %dx -> assembly likely fails", c, mg);
      if(n+0 < gn) r = r sprintf("; read N50 %s bp < %d bp (short)", n, gn);
      print r }')

  # ---- accumulate summary row + report arrays ----
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$n_all" "$bp_all" "$n50" "$gc" "$pct_its" "$cov_all" "$cov_1k" "$verdict" "$top_species" \
    >> "$SUMMARY"

  R_SAMPLE+=("$SAMPLE");  R_NALL+=("$n_all");   R_BPALL+=("$bp_all")
  R_N50+=("$n50");        R_GC+=("$gc");        R_PITS+=("$pct_its")
  R_COVALL+=("$cov_all"); R_COV1K+=("$cov_1k"); R_VERDICT+=("$verdict")
  R_AMP+=("$amp_class");  R_TOP+=("$top_species"); R_REASON+=("$reason")
  R_B+=("$(printf '%s' "$BUCKETS" | tr '\t' ',')")
done

###############################################################################
# CHECK 8 — ARTIFACTS
###############################################################################
banner "CHECK 8  ARTIFACTS"

echo "Per-sample summary:"
column -t -s "$(printf '\t')" "$SUMMARY"
echo

# fungiforge samplesheet stub
STUB="$OUTDIR_ABS/samplesheet.stub.csv"
{
  echo "sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season"
  for f in "${SAMPLE_FILES[@]}"; do
    s="$(sample_of "$f")"
    printf '%s,%s,NA,NA,NA,NA,NA\n' "$s" "$(abspath "$f")"
  done
} > "$STUB"

###############################################################################
# REPORT — saveable Markdown + self-contained HTML, with figures
###############################################################################
sanitize(){ printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'; }

# ---- run-overview numbers ----
NSAMP=${#R_SAMPLE[@]}
TOTAL_BP=0; for v in "${R_BPALL[@]}"; do TOTAL_BP=$(( TOTAL_BP + ${v:-0} )); done
TOTAL_H=$(awk -v b="$TOTAL_BP" 'BEGIN{ if(b>=1e9)printf "%.2f Gb",b/1e9; else if(b>=1e6)printf "%.1f Mb",b/1e6; else printf "%d bp",b}')

# ---- dump the data the plotter consumes ----
DATATSV="$OUTDIR_ABS/.scripts/report_data.tsv"
{ printf 'sample\tn_all\tbp_all\tn50\tgc\tpct_its\tcov_all\tcov1k\tverdict\tamp\ttop\tbuckets\n'
  for i in "${!R_SAMPLE[@]}"; do
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${R_SAMPLE[$i]}" "${R_NALL[$i]}" "${R_BPALL[$i]}" "${R_N50[$i]}" "${R_GC[$i]}" \
      "${R_PITS[$i]}" "${R_COVALL[$i]}" "${R_COV1K[$i]}" "${R_VERDICT[$i]}" "${R_AMP[$i]}" \
      "${R_TOP[$i]}" "${R_B[$i]}"
  done
} > "$DATATSV"

# ---- figures (matplotlib in-container; degrade to ASCII if unavailable) ----
FIGDIR="$OUTDIR_ABS/figs"; mkdir -p "$FIGDIR"
cat > "$OUTDIR_ABS/.scripts/mkfigs.py" <<'PY'
import os, sys, csv, re
try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as e:
    sys.stderr.write("no-matplotlib: %s\n" % e); sys.exit(3)

DATA="/out/.scripts/report_data.tsv"; FIG="/out/figs"; os.makedirs(FIG, exist_ok=True)
COV_GO=float(os.environ.get("COV_GO","30"))
GOOD_N50=float(os.environ.get("GOOD_N50","5000"))
AMP_HI=float(os.environ.get("AMP_PCT_HI","50"))
def safe(s): return re.sub(r"[^A-Za-z0-9_.-]","_",s)

rows=[]
with open(DATA) as fh:
    for r in csv.DictReader(fh, delimiter="\t"): rows.append(r)
if not rows: sys.exit(0)
samples=[x["sample"] for x in rows]; idx=np.arange(len(samples))
def num(x):
    try: return float(x)
    except: return 0.0
def save(fig,name): fig.savefig(os.path.join(FIG,name),dpi=120,bbox_inches="tight"); plt.close(fig); print(name)

labels=["<300","300-500","500-800","800-1200","1200-3000",">3000"]
for x in rows:
    b=[int(v) for v in x["buckets"].split(",")] if x.get("buckets") else [0]*6
    fig,ax=plt.subplots(figsize=(6,3.2))
    ax.bar(range(6),b,color="#4C78A8")
    ax.set_xticks(range(6)); ax.set_xticklabels(labels,rotation=30,ha="right",fontsize=8)
    ax.set_ylabel("reads"); ax.set_title("Read-length distribution — %s"%x["sample"])
    save(fig,"hist_%s.png"%safe(x["sample"]))

cov_all=[num(x["cov_all"]) for x in rows]; cov1k=[num(x["cov1k"]) for x in rows]
fig,ax=plt.subplots(figsize=(max(6,len(samples)*0.9),3.6)); w=0.38
ax.bar(idx-w/2,cov_all,w,label="all reads",color="#9ecae1")
ax.bar(idx+w/2,cov1k,w,label=">=1kb reads",color="#3182bd")
ax.axhline(COV_GO,color="#e6550d",ls="--",lw=1,label="GO threshold %gx"%COV_GO)
ax.set_xticks(idx); ax.set_xticklabels(samples,rotation=30,ha="right",fontsize=8)
ax.set_ylabel("fold coverage"); ax.set_title("Estimated coverage per sample"); ax.legend(fontsize=8)
save(fig,"coverage.png")

n50=[num(x["n50"]) for x in rows]
fig,ax=plt.subplots(figsize=(max(6,len(samples)*0.9),3.6))
ax.bar(idx,n50,color="#54a24b")
ax.axhline(GOOD_N50,color="#e6550d",ls="--",lw=1,label="good ONT %g bp"%GOOD_N50)
ax.set_xticks(idx); ax.set_xticklabels(samples,rotation=30,ha="right",fontsize=8)
ax.set_ylabel("read N50 (bp)"); ax.set_title("Read N50 per sample"); ax.legend(fontsize=8)
save(fig,"read_n50.png")

pits=[num(x["pct_its"]) for x in rows]
fig,ax=plt.subplots(figsize=(max(6,len(samples)*0.9),3.6))
ax.bar(idx,pits,color="#b279a2")
ax.axhline(AMP_HI,color="#e6550d",ls="--",lw=1,label="amplicon flag %g%%"%AMP_HI)
ax.set_xticks(idx); ax.set_xticklabels(samples,rotation=30,ha="right",fontsize=8)
ax.set_ylabel("% subsampled reads hitting ITS"); ax.set_title("ITS-hit fraction (amplicon vs WGS)"); ax.legend(fontsize=8)
save(fig,"pct_its.png")
PY

export COV_GO GOOD_N50 AMP_PCT_HI
echo "[..] rendering figures (matplotlib)"
FIG_LOG="$(dock '/opt/conda/bin/python /out/.scripts/mkfigs.py 2>&1' || true)"
shopt -s nullglob; _figs=("$FIGDIR"/*.png); shopt -u nullglob
if [ "${#_figs[@]}" -gt 0 ]; then
  FIGS_OK=1; echo "[ok] wrote ${#_figs[@]} figure(s) to $FIGDIR"
else
  FIGS_OK=0; warn "figure rendering unavailable — reports will use ASCII bars. Detail: $(printf '%s' "$FIG_LOG" | tail -1)"
fi

# ---- helpers for report bodies ----
ascii_bar(){ awk -v v="$1" -v m="$2" -v w="${3:-38}" 'BEGIN{ if(m<=0)m=1; n=int(w*v/m+0.5); if(n>w)n=w; if(n<0)n=0; s=""; for(i=0;i<n;i++)s=s"#"; printf "%s", s }'; }
maxof(){ awk 'BEGIN{m=0} {if($1+0>m)m=$1} END{print (m>0?m:1)}'; }
b64img(){ printf 'data:image/png;base64,'; base64 < "$1" | tr -d '\n'; }

# ---- category lists for recommendations ----
GO_L=(); MARG_L=(); NOGO_L=(); AMP_L=(); RESEQ_L=()
for i in "${!R_SAMPLE[@]}"; do
  s="${R_SAMPLE[$i]}"
  case "${R_VERDICT[$i]}" in
    GO)        GO_L+=("$s") ;;
    MARGINAL)  MARG_L+=("$s") ;;
    NO-GO)     NOGO_L+=("$s") ;;
  esac
  [ "${R_AMP[$i]}" = "AMPLICON" ] && AMP_L+=("$s")
  if [ "${R_VERDICT[$i]}" = "NO-GO" ] || awk -v n="${R_N50[$i]}" -v g="$GOOD_N50" 'BEGIN{exit !(n+0<g)}'; then
    RESEQ_L+=("$s")
  fi
done
# comma-join, tolerant of empty arrays under `set -u` (macOS bash 3.2 safe)
listjoin(){ eval "local c=\${#$1[@]}"; if [ "$c" -eq 0 ]; then printf '(none)'; return; fi
  eval "local arr=(\"\${$1[@]}\")"; local out="" e
  for e in "${arr[@]}"; do [ -n "$out" ] && out="$out, $e" || out="$e"; done; printf '%s' "$out"; }
GO_S="$(listjoin GO_L)"; MARG_S="$(listjoin MARG_L)"; NOGO_S="$(listjoin NOGO_L)"
AMP_S="$(listjoin AMP_L)"; RESEQ_S="$(listjoin RESEQ_L)"

# maxima for ASCII bars
COVMAX=$(printf '%s\n' "${R_COVALL[@]}" "${R_COV1K[@]}" "$COV_GO" | maxof)
N50MAX=$(printf '%s\n' "${R_N50[@]}" "$GOOD_N50" | maxof)

###############################################################################
# Markdown report
###############################################################################
MDREPORT="$OUTDIR_ABS/preflight_qc_report.md"
{
  echo "# FungiForge preflight QC report"
  echo
  echo "**Generated:** $RUN_DATE on \`$RUN_HOST\`"
  echo
  echo "## Run overview"
  echo
  echo "| field | value |"
  echo "|---|---|"
  echo "| input path | \`$INPUT\` |"
  echo "| database | \`$DB\` (UNITE: \`$UNITE\`) |"
  echo "| base image | \`$IMAGE\` |"
  echo "| merged reads | \`$MERGED_ABS\` |"
  echo "| samples | $NSAMP |"
  echo "| total data | $TOTAL_H |"
  echo "| genome size (coverage) | $GENOME_SIZE bp |"
  echo
  echo "## Environment readiness"
  echo
  echo "| check | status | detail |"
  echo "|---|---|---|"
  for row in "${ENV_ROWS[@]}"; do
    IFS='|' read -r st lb dt <<<"$row"
    echo "| $lb | **$st** | $dt |"
  done
  echo
  echo "## Per-sample QC summary"
  echo
  echo "| sample | reads | sum_bp | read_N50 | GC% | %ITS | cov_all | cov_ge1kb | verdict | top species (read-based ID) |"
  echo "|---|--:|--:|--:|--:|--:|--:|--:|:--:|---|"
  for i in "${!R_SAMPLE[@]}"; do
    echo "| ${R_SAMPLE[$i]} | ${R_NALL[$i]} | ${R_BPALL[$i]} | ${R_N50[$i]} | ${R_GC[$i]} | ${R_PITS[$i]} | ${R_COVALL[$i]}x | ${R_COV1K[$i]}x | ${R_VERDICT[$i]} | ${R_TOP[$i]} |"
  done
  echo
  echo "## Figures"
  echo
  if [ "$FIGS_OK" -eq 1 ]; then
    echo "![Coverage per sample](figs/coverage.png)"
    echo
    echo "![Read N50 per sample](figs/read_n50.png)"
    echo
    echo "![ITS-hit fraction](figs/pct_its.png)"
    echo
    echo "**Per-sample read-length distributions:**"
    echo
    for i in "${!R_SAMPLE[@]}"; do
      sf="$(sanitize "${R_SAMPLE[$i]}")"
      echo "![${R_SAMPLE[$i]} read lengths](figs/hist_${sf}.png)"
      echo
    done
  else
    echo "_matplotlib was unavailable; ASCII bar charts shown instead._"
    echo
    echo '```'
    echo "Coverage >=1kb (bar) vs GO threshold ${COV_GO}x:"
    for i in "${!R_SAMPLE[@]}"; do
      printf '  %-12s %6sx |%s\n' "${R_SAMPLE[$i]}" "${R_COV1K[$i]}" "$(ascii_bar "${R_COV1K[$i]}" "$COVMAX")"
    done
    echo
    echo "Read N50 (bar) vs good-ONT ${GOOD_N50} bp:"
    for i in "${!R_SAMPLE[@]}"; do
      printf '  %-12s %8s |%s\n' "${R_SAMPLE[$i]}" "${R_N50[$i]}" "$(ascii_bar "${R_N50[$i]}" "$N50MAX")"
    done
    echo '```'
  fi
  echo
  echo "## Amplicon-vs-WGS determination"
  echo
  for i in "${!R_SAMPLE[@]}"; do
    case "${R_AMP[$i]}" in
      AMPLICON)  msg="looks like an **ITS amplicon** library (${R_PITS[$i]}% of subsampled reads hit ITS) — fungiforge is the WRONG tool; route to an ITS/metabarcoding workflow." ;;
      WGS)       msg="looks like **whole-genome shotgun** (${R_PITS[$i]}% ITS hits) — correct input for fungiforge." ;;
      *)         msg="**ambiguous** (${R_PITS[$i]}% ITS hits) — inspect the read-length distribution and species spread." ;;
    esac
    echo "- **${R_SAMPLE[$i]}**: $msg"
  done
  echo
  echo "## Recommendations"
  echo
  echo "Verdict thresholds (heuristic, on \`>=1kb\` fold-coverage): **GO** \`>= ${COV_GO}x\`, **MARGINAL** \`${COV_MARGINAL}-${COV_GO}x\`, **NO-GO** \`< ${COV_MARGINAL}x\`. A read N50 below **${GOOD_N50} bp** is flagged as short."
  echo
  echo "### Per sample"
  echo
  for i in "${!R_SAMPLE[@]}"; do
    echo "- **${R_SAMPLE[$i]} — ${R_VERDICT[$i]}**: ${R_REASON[$i]}. Top read-based ID: ${R_TOP[$i]}."
  done
  echo
  echo "### How to proceed"
  echo
  echo "- **Attempt de novo assembly:** ${GO_S}."
  echo "- **Marginal — assemble only with relaxed settings** (e.g. lower \`--ont_min_len\` to admit more reads; expect a fragmented, draft-quality assembly): ${MARG_S}."
  echo "- **ID-only (assembly not worthwhile):** ${NOGO_S}. Use the read-based species ID above; do not spend the multi-hour run hoping for a genome."
  echo "- **Consider re-sequencing** (short read-N50 and/or low yield): ${RESEQ_S}."
  echo "- **Not shotgun — do NOT run fungiforge** (amplicon libraries): ${AMP_S}."
  echo
  echo "---"
  echo
  echo "_Verdicts are heuristics for triage, not guarantees. Yeasts (~15 Mb) have roughly double the coverage implied here (default genome size assumes a ~${GENOME_SIZE} bp mold). ONT resistance/point-mutation calls remain provisional until hybrid-polished._"
} > "$MDREPORT"

###############################################################################
# Self-contained HTML report (figures embedded as base64 data-URIs)
###############################################################################
HTMLREPORT="$OUTDIR_ABS/preflight_qc_report.html"
{
  cat <<'HEAD'
<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>FungiForge preflight QC report</title>
<style>
:root{color-scheme:light dark}
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;line-height:1.5;
     max-width:1000px;margin:2rem auto;padding:0 1rem;color:#1a1a1a;background:#fff}
h1,h2,h3{line-height:1.25} h1{border-bottom:2px solid #3182bd;padding-bottom:.3rem}
h2{margin-top:2rem;border-bottom:1px solid #ccc;padding-bottom:.2rem}
table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.92rem;display:block;overflow-x:auto}
th,td{border:1px solid #ccc;padding:.35rem .55rem;text-align:left}
th{background:#eef4fa} code{background:#f0f0f0;padding:.05rem .3rem;border-radius:3px}
img{max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px;margin:.4rem 0}
.GO,.PASS{color:#0a7d29;font-weight:700}.MARGINAL,.WARN{color:#b26a00;font-weight:700}
.NO-GO,.FAIL{color:#c0261a;font-weight:700}
.rec{background:#f6fafd;border-left:4px solid #3182bd;padding:.6rem 1rem;border-radius:4px}
pre{background:#f0f0f0;padding:.7rem;border-radius:4px;overflow-x:auto}
@media(prefers-color-scheme:dark){
 body{color:#e6e6e6;background:#161616}th{background:#1e2a36}code,pre{background:#222}
 th,td{border-color:#444}img{border-color:#444}.rec{background:#12202c}}
</style></head><body>
HEAD
  echo "<h1>FungiForge preflight QC report</h1>"
  echo "<p><strong>Generated:</strong> $RUN_DATE on <code>$RUN_HOST</code></p>"

  echo "<h2>Run overview</h2><table>"
  echo "<tr><th>input path</th><td><code>$INPUT</code></td></tr>"
  echo "<tr><th>database</th><td><code>$DB</code> (UNITE: <code>$UNITE</code>)</td></tr>"
  echo "<tr><th>base image</th><td><code>$IMAGE</code></td></tr>"
  echo "<tr><th>merged reads</th><td><code>$MERGED_ABS</code></td></tr>"
  echo "<tr><th>samples</th><td>$NSAMP</td></tr>"
  echo "<tr><th>total data</th><td>$TOTAL_H</td></tr>"
  echo "<tr><th>genome size (coverage)</th><td>$GENOME_SIZE bp</td></tr></table>"

  echo "<h2>Environment readiness</h2><table><tr><th>check</th><th>status</th><th>detail</th></tr>"
  for row in "${ENV_ROWS[@]}"; do
    IFS='|' read -r st lb dt <<<"$row"
    echo "<tr><td>$lb</td><td class=\"$st\">$st</td><td>$dt</td></tr>"
  done
  echo "</table>"

  echo "<h2>Per-sample QC summary</h2><table>"
  echo "<tr><th>sample</th><th>reads</th><th>sum_bp</th><th>read_N50</th><th>GC%</th><th>%ITS</th><th>cov_all</th><th>cov_ge1kb</th><th>verdict</th><th>top species</th></tr>"
  for i in "${!R_SAMPLE[@]}"; do
    echo "<tr><td>${R_SAMPLE[$i]}</td><td>${R_NALL[$i]}</td><td>${R_BPALL[$i]}</td><td>${R_N50[$i]}</td><td>${R_GC[$i]}</td><td>${R_PITS[$i]}</td><td>${R_COVALL[$i]}x</td><td>${R_COV1K[$i]}x</td><td class=\"${R_VERDICT[$i]}\">${R_VERDICT[$i]}</td><td>${R_TOP[$i]}</td></tr>"
  done
  echo "</table>"

  echo "<h2>Figures</h2>"
  if [ "$FIGS_OK" -eq 1 ]; then
    for fig in coverage read_n50 pct_its; do
      [ -f "$FIGDIR/$fig.png" ] && echo "<img alt=\"$fig\" src=\"$(b64img "$FIGDIR/$fig.png")\">"
    done
    echo "<h3>Per-sample read-length distributions</h3>"
    for i in "${!R_SAMPLE[@]}"; do
      sf="$(sanitize "${R_SAMPLE[$i]}")"
      [ -f "$FIGDIR/hist_${sf}.png" ] && echo "<img alt=\"${R_SAMPLE[$i]} read lengths\" src=\"$(b64img "$FIGDIR/hist_${sf}.png")\">"
    done
  else
    echo "<p><em>matplotlib unavailable — ASCII bars:</em></p><pre>"
    echo "Coverage >=1kb vs GO threshold ${COV_GO}x:"
    for i in "${!R_SAMPLE[@]}"; do
      printf '  %-12s %6sx |%s\n' "${R_SAMPLE[$i]}" "${R_COV1K[$i]}" "$(ascii_bar "${R_COV1K[$i]}" "$COVMAX")"
    done
    echo
    echo "Read N50 vs good-ONT ${GOOD_N50} bp:"
    for i in "${!R_SAMPLE[@]}"; do
      printf '  %-12s %8s |%s\n' "${R_SAMPLE[$i]}" "${R_N50[$i]}" "$(ascii_bar "${R_N50[$i]}" "$N50MAX")"
    done
    echo "</pre>"
  fi

  echo "<h2>Amplicon-vs-WGS determination</h2><ul>"
  for i in "${!R_SAMPLE[@]}"; do
    case "${R_AMP[$i]}" in
      AMPLICON)  msg="looks like an <strong>ITS amplicon</strong> library (${R_PITS[$i]}% ITS hits) &mdash; fungiforge is the WRONG tool; route to an ITS/metabarcoding workflow." ;;
      WGS)       msg="looks like <strong>whole-genome shotgun</strong> (${R_PITS[$i]}% ITS hits) &mdash; correct input for fungiforge." ;;
      *)         msg="<strong>ambiguous</strong> (${R_PITS[$i]}% ITS hits) &mdash; inspect read lengths and species spread." ;;
    esac
    echo "<li><strong>${R_SAMPLE[$i]}</strong>: $msg</li>"
  done
  echo "</ul>"

  echo "<h2>Recommendations</h2>"
  echo "<p>Verdict thresholds (heuristic, on <code>&ge;1kb</code> fold-coverage): <span class=\"GO\">GO</span> &ge; ${COV_GO}x, <span class=\"MARGINAL\">MARGINAL</span> ${COV_MARGINAL}&ndash;${COV_GO}x, <span class=\"NO-GO\">NO-GO</span> &lt; ${COV_MARGINAL}x. Read N50 below ${GOOD_N50} bp is flagged short.</p>"
  echo "<h3>Per sample</h3><ul>"
  for i in "${!R_SAMPLE[@]}"; do
    echo "<li><strong>${R_SAMPLE[$i]} &mdash; <span class=\"${R_VERDICT[$i]}\">${R_VERDICT[$i]}</span></strong>: ${R_REASON[$i]}. Top read-based ID: ${R_TOP[$i]}.</li>"
  done
  echo "</ul><div class=\"rec\"><h3>How to proceed</h3><ul>"
  echo "<li><strong>Attempt de novo assembly:</strong> ${GO_S}.</li>"
  echo "<li><strong>Marginal &mdash; assemble only with relaxed settings</strong> (e.g. lower <code>--ont_min_len</code>; expect a fragmented draft): ${MARG_S}.</li>"
  echo "<li><strong>ID-only (assembly not worthwhile):</strong> ${NOGO_S}. Use the read-based ID; don't spend the multi-hour run on these.</li>"
  echo "<li><strong>Consider re-sequencing</strong> (short read-N50 / low yield): ${RESEQ_S}.</li>"
  echo "<li><strong>Not shotgun &mdash; do NOT run fungiforge</strong> (amplicon): ${AMP_S}.</li>"
  echo "</ul></div>"
  echo "<p><em>Verdicts are heuristics for triage, not guarantees. Yeasts (~15 Mb) have ~double the coverage implied here. ONT point-mutation calls remain provisional until hybrid-polished.</em></p>"
  echo "</body></html>"
} > "$HTMLREPORT"

# tidy the staged helper dir (keep figs/, drop scratch)
rm -rf "$OUTDIR_ABS/.scripts" 2>/dev/null || true

echo "Wrote:"
echo "  summary TSV      : $SUMMARY"
echo "  QC report (md)   : $MDREPORT"
echo "  QC report (html) : $HTMLREPORT   (self-contained; figures embedded)"
[ "$FIGS_OK" -eq 1 ] && echo "  figures (png)    : $FIGDIR/"
echo "  samplesheet stub : $STUB   (fill compartment/facility/season, add Illumina if hybrid)"
echo "  merged reads     : $MERGED_ABS"
echo
echo "Next: edit the stub, then run fungiforge, e.g."
echo "  nextflow run main.nf -profile local,docker --samplesheet '$STUB' --data_dir '$DB'"
echo
echo "Reminder: only samples with a GO/MARGINAL verdict (check 6) are worth assembling;"
echo "NO-GO samples still get an identity from check 7 (read-based ID)."
