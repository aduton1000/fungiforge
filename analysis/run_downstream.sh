#!/usr/bin/env bash
# fungiforge comparative layer — run the objective scripts in order against the
# per-isolate master tables produced by Layer 1 (mirrors the bacterial run_downstream.sh).
#
#   analysis/run_downstream.sh [RESULTS_DIR] [SAMPLE_METADATA_CSV]
#
# RESULTS_DIR defaults to ../results ; SAMPLE_METADATA_CSV is optional (air volume,
# health-seeking, etc.) joined by sample id. Requires Rscript + the fungiforge R deps
# (dplyr, readr, tidyr, ggplot2, stringr, vegan, lme4, broom.mixed).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RS="${RSCRIPT:-Rscript}"
RESULTS="${1:-$HERE/../results}"
META="${2:-}"
LOG="$HERE/output/downstream_$(date +%Y%m%d_%H%M%S).log"; mkdir -p "$HERE/output"

run(){ echo ">>> $*" | tee -a "$LOG"; "$RS" "$@" 2>&1 | tee -a "$LOG"; }

echo "== fungiforge downstream · results=$RESULTS ==" | tee "$LOG"
run "$HERE/scripts/01_merge_metadata.R" "$RESULTS" ${META:+"$META"}
run "$HERE/scripts/02_recover_metadata.R"
for s in 03_objective1_composition 04_objective2_resistome 05_objective3_mobile \
         06_objective4_transmission 07_objective5_risk 08_objective6_virulence \
         09_objective7_connectivity 10_objective8_novelty; do
  run "$HERE/scripts/$s.R"
done
echo "== done · tables in output/ · figures in figures/ ==" | tee -a "$LOG"
