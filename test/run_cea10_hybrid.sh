#!/usr/bin/env bash
# CEA10 A. fumigatus HYBRID end-to-end test (ONT SRR28036069 + Illumina SRR13127874)
set -euo pipefail
REPO="~/fungiforge"
DB="/Volumes/Extreme SSD/fungiforge_db"
WORK="/Volumes/Extreme SSD/fungiforge_work_cea10"
cd "$REPO"
# stale-lock hygiene
pkill -9 -f nextflow.cli.Launcher 2>/dev/null || true
find .nextflow -name LOCK -o -name '*.lock' 2>/dev/null | xargs -r rm -f
BGC_FLAG="--skip_bgc"
[ -f "$DB/antismash/.done" ] && BGC_FLAG=""   # enable BGC only when antiSMASH DB is ready
echo "[run] BGC flag: '${BGC_FLAG:-<enabled>}'"
exec conda run --no-capture-output -n nextflow nextflow run main.nf \
  -profile local,docker \
  --samplesheet test/sra/samplesheet.AFUM_CEA10_HYBRID.csv \
  --data_dir "$DB" --outdir results_cea10 \
  --busco_lineage fungi_odb10 \
  $BGC_FLAG \
  -w "$WORK" -resume
