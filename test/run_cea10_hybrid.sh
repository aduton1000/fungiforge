#!/usr/bin/env bash
# Example: A. fumigatus CEA10 hybrid end-to-end test (ONT SRR28036069 + Illumina SRR13127874).
# Fetch the reads with test/fetch_sra_test.sh, drop a samplesheet at the path below, then run.
#
# IMPORTANT (portability lessons baked in):
#   * WORK, FUNANNOTATE_DB and ANTISMASH_DB MUST be SPACE-FREE paths. Flye aborts on spaces in
#     read paths; funannotate's embedded BUSCO2 and antiSMASH pass DB paths to Augustus/blastp
#     UNQUOTED, so a space (e.g. an external drive named with a space) silently yields 0 genes /
#     0 clusters. Reference DBs may live on a space-containing drive, but point these three at a
#     space-free location (copy funannotate/ + antismash/ there if your DB root has spaces).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
DB="${FUNGIFORGE_DB:?set FUNGIFORGE_DB to your reference-database root}"
WORK="${FUNGIFORGE_WORK:-$HOME/fungiforge_work_cea10}"          # space-free
FUNDB="${FUNANNOTATE_DB_DIR:-$HOME/fungiforge_fundb}"           # space-free copy of $DB/funannotate
ASDB="${ANTISMASH_DB_DIR:-$HOME/fungiforge_asdb}"               # space-free copy of $DB/antismash
SHEET="${SAMPLESHEET:-$REPO/test/sra/samplesheet.AFUM_CEA10_HYBRID.csv}"
cd "$REPO"
# stale-lock hygiene
pkill -9 -f nextflow.cli.Launcher 2>/dev/null || true
find .nextflow -name LOCK -o -name '*.lock' 2>/dev/null | xargs -r rm -f
BGC_FLAG="--skip_bgc"
[ -f "$ASDB/.done" ] || [ -f "$DB/antismash/.done" ] && BGC_FLAG=""   # enable BGC once the antiSMASH DB is prebuilt
echo "[run] BGC flag: '${BGC_FLAG:-<enabled>}'"
exec conda run --no-capture-output -n nextflow nextflow run main.nf \
  -profile local,docker \
  --samplesheet "$SHEET" \
  --data_dir "$DB" --outdir results_cea10 \
  --busco_lineage fungi_odb10 \
  --funannotate_db "$FUNDB" \
  --antismash_db "$ASDB" \
  $BGC_FLAG \
  -w "$WORK" -resume
