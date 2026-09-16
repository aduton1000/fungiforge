#!/usr/bin/env bash
# Fetch the public control isolates used by the validation suite (W1.1) into
#   $FUNGIFORGE_DB/controls/<name>/   (reads + a samplesheet row; `.done` marker per control)
# and verify each control's genotype FROM ITS READS before it is trusted (bin/control_genotype.py
# against test/controls/*.json). A control whose reads do not carry the claimed genotype is
# marked FAILED_GENOTYPE and must not be used — this is what the deposited "C87" reads taught us.
#
#   bash bin/fetch_controls.sh [control ...]        default: all controls below
#
# Controls (A. fumigatus cyp51A; ONT, PromethION, Dorado hac 400 bps; PRJNA1494181, 2026):
#   afum_TR34_L98H_226BB2   SRR39591703  TR34 + L98H (+ F46Y M172V N284T D255E E427K background)
#   afum_TR34_L98H_160CME3  SRR39591704  TR34 + L98H   (second positive)
#   afum_wt_157DB3          SRR39591705  no TR, no L98H, same F46Y/M172V/N284T/D255E/E427K background
#   afum_CEA10              SRR28036069 (ONT) + SRR13127874 (Illumina)  wild-type reference isolate
# Downloads use ENA's FASTQ mirror when the run is mirrored, else `fasterq-dump` (sra-tools)
# if installed, else NCBI's FASTQ streaming endpoint (no tools needed, slower). Reads are 4–8 GB
# per ONT run. The read-level check needs minimap2: set FUNGIFORGE_TOOL_RUNNER to a command
# prefix that provides it when it is not on PATH, e.g. on the cluster
#   FUNGIFORGE_TOOL_RUNNER="apptainer exec -B /hpc,$HOME /hpc/opt/fungiforge-dev/images/fungiforge-0.1.0.sif"
set -uo pipefail
DB="${FUNGIFORGE_DB:?set FUNGIFORGE_DB to the reference-database root}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CDIR="$DB/controls"; LOGDIR="$DB/logs"; mkdir -p "$CDIR" "$LOGDIR"
MANIFEST="$CDIR/MANIFEST.tsv"; [ -f "$MANIFEST" ] || echo -e "control\truns\tgenotype_check\tstatus\ttimestamp" > "$MANIFEST"
log(){ echo "[$(date '+%F %T')] $*"; }
note(){ echo -e "$1\t$2\t$3\t$4\t$(date '+%F %T')" >> "$MANIFEST"; }

dl(){  # dl <url> <out>
  if command -v aria2c >/dev/null 2>&1; then
    aria2c --continue=true -x8 -s8 -k1M --max-tries=0 --retry-wait=10 --file-allocation=none \
           --auto-file-renaming=false --console-log-level=warn -o "$2" "$1"
  else
    curl -L --retry 10 --retry-delay 10 -C - -o "$2" "$1"
  fi
}

# ENA mirror URLs for a run (may be empty when not yet mirrored)
ena_fastq_urls(){
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp" \
    | awk -F'\t' 'NR>1 && $2!="" {print $2}' | tr ';' '\n' | sed 's#^#https://#'
}

fetch_run(){  # fetch_run <run> <dir>  -> files <dir>/<run>*.fastq.gz
  local run="$1" dir="$2" urls
  mkdir -p "$dir"
  if ls "$dir/$run"*.fastq.gz >/dev/null 2>&1; then log "have $run — skip"; return 0; fi
  urls="$(ena_fastq_urls "$run")"
  if [ -n "$urls" ]; then
    local u; for u in $urls; do ( cd "$dir" && dl "$u" "$(basename "$u")" ) || return 1; done
  elif command -v fasterq-dump >/dev/null 2>&1; then
    log "$run not on the ENA mirror — fasterq-dump"
    ( cd "$dir" && fasterq-dump --split-files --threads 8 -O . "$run" && gzip -f "$run"*.fastq ) || return 1
  else
    log "$run not on the ENA mirror and sra-tools not installed — streaming FASTQ from NCBI (slow)"
    ( cd "$dir" && curl -L --retry 10 --retry-delay 30 -sS "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=$run" \
        | gzip > "$run.fastq.gz.part" && [ -s "$run.fastq.gz.part" ] && gzip -t "$run.fastq.gz.part" && mv "$run.fastq.gz.part" "$run.fastq.gz" ) \
      || { rm -f "$dir/$run.fastq.gz.part"; return 1; }
  fi
}

verify(){  # verify <control> <spec-control-name> <platform> <reads...>
  local name="$1" ctrl="$2" platform="$3"; shift 3
  local args=(); local r; for r in "$@"; do args+=(--reads "$r"); done
  # shellcheck disable=SC2086  # the runner prefix is a command with arguments
  ${FUNGIFORGE_TOOL_RUNNER:-} python3 "$REPO_DIR/bin/control_genotype.py" --spec "$REPO_DIR/test/controls/afum_cyp51A_controls.json" \
      --control "$ctrl" --platform "$platform" --threads 8 --out "$CDIR/$name/genotype_check.json" "${args[@]}" \
      >"$CDIR/$name/genotype_check.log" 2>&1
}

sheet_row(){  # sheet_row <control> <sample> <ont> <r1> <r2>
  printf 'sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season\n%s,%s,%s,%s,control,public,NA\n' \
    "$2" "$3" "$4" "$5" > "$CDIR/$1/samples.csv"
}

control_ont(){  # control_ont <name> <run> <spec-control>
  local name="$1" run="$2" ctrl="$3" d="$CDIR/$1"
  [ -f "$d/.done" ] && { log "$name done — skip"; return; }
  fetch_run "$run" "$d" || { note "$name" "$run" "-" "FAILED_DOWNLOAD"; return; }
  local reads; reads="$(ls "$d/$run"*.fastq.gz | head -1)"
  if verify "$name" "$ctrl" ont "$reads"; then
    sheet_row "$name" "$name" "$reads" "" ""; date '+%F %T' > "$d/.done"; note "$name" "$run" "PASS" "OK"; log "✓ $name ($run): genotype verified from reads"
  else
    note "$name" "$run" "FAIL" "FAILED_GENOTYPE"; log "✗ $name ($run): reads do NOT carry the expected genotype — see $d/genotype_check.log"
  fi
}

control_hybrid(){  # control_hybrid <name> <ont-run> <illumina-run> <spec-control>
  local name="$1" ont="$2" ilmn="$3" ctrl="$4" d="$CDIR/$1"
  [ -f "$d/.done" ] && { log "$name done — skip"; return; }
  fetch_run "$ont" "$d" && fetch_run "$ilmn" "$d" || { note "$name" "$ont;$ilmn" "-" "FAILED_DOWNLOAD"; return; }
  local o r1 r2; o="$(ls "$d/$ont"*.fastq.gz | head -1)"; r1="$(ls "$d/$ilmn"_1.fastq.gz)"; r2="$(ls "$d/$ilmn"_2.fastq.gz)"
  if verify "$name" "$ctrl" ont "$o" && verify "${name}" "$ctrl" illumina "$r1" "$r2"; then
    sheet_row "$name" "$name" "$o" "$r1" "$r2"; date '+%F %T' > "$d/.done"; note "$name" "$ont;$ilmn" "PASS" "OK"; log "✓ $name: genotype verified from both read sets"
  else
    note "$name" "$ont;$ilmn" "FAIL" "FAILED_GENOTYPE"; log "✗ $name: reads do NOT carry the expected genotype — see $d/genotype_check.log"
  fi
}

run_control(){
  case "$1" in
    afum_TR34_L98H_226BB2)  control_ont "$1" SRR39591703 TR34_L98H ;;
    afum_TR34_L98H_160CME3) control_ont "$1" SRR39591704 TR34_L98H ;;
    afum_wt_157DB3)         control_ont "$1" SRR39591705 wild_type ;;
    afum_CEA10)             control_hybrid "$1" SRR28036069 SRR13127874 wild_type ;;
    *) log "unknown control: $1"; return 1 ;;
  esac
}

CONTROLS=("$@"); [ ${#CONTROLS[@]} -eq 0 ] && CONTROLS=(afum_TR34_L98H_226BB2 afum_wt_157DB3 afum_TR34_L98H_160CME3 afum_CEA10)
log "==== fetch_controls · $CDIR · ${CONTROLS[*]} ===="
for c in "${CONTROLS[@]}"; do run_control "$c"; done
log "manifest: $MANIFEST"; column -t -s $'\t' "$MANIFEST" 2>/dev/null || cat "$MANIFEST"
