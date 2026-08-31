#!/usr/bin/env bash
# ==============================================================================
# fungiforge — fetch a real public fungal dataset (ONT + Illumina) from ENA and
# build a samplesheet, to test the pipeline end-to-end on genuine data.
#
#   bin/fetch_sra_test.sh --sample TEST_SCER --ont ERRxxxx --r1r2 ERRyyyy \
#        --compartment AIR --outdir test/sra
#
# Resolves ENA fastq URLs from run accessions (filereport API), downloads with
# aria2c (resumable), optionally subsamples ONT for a fast pass, and writes
# samplesheet.csv ready for `nextflow run main.nf --samplesheet ...`.
#
# Candidate datasets to run (verify current accessions at run time via ENA search):
#   * Saccharomyces cerevisiae (~12 Mb) — fast full pass / CI
#   * Aspergillus fumigatus — exercises the cyp51A TR34/TR46 + azole panel
#   * Candida auris — FKS/ERG11 + clade ID
# ==============================================================================
set -uo pipefail
SAMPLE=""; ONT=""; PAIRED=""; COMP="AIR"; FAC="TEST"; SEASON="Dry"; OUT="test/sra"; SUBSAMPLE=""
while [ $# -gt 0 ]; do case "$1" in
  --sample) SAMPLE="$2"; shift 2;;
  --ont)    ONT="$2"; shift 2;;             # ONT run accession (ERR/SRR/DRR)
  --r1r2)   PAIRED="$2"; shift 2;;          # Illumina paired run accession
  --compartment) COMP="$2"; shift 2;;
  --facility) FAC="$2"; shift 2;;
  --season) SEASON="$2"; shift 2;;
  --outdir) OUT="$2"; shift 2;;
  --subsample) SUBSAMPLE="$2"; shift 2;;    # e.g. 200000 ONT reads for a fast pass
  *) echo "unknown arg $1"; exit 1;;
esac; done
[ -z "$SAMPLE" ] && { echo "need --sample"; exit 1; }
mkdir -p "$OUT/$SAMPLE"

dl(){ if command -v aria2c >/dev/null 2>&1; then
        aria2c --continue=true -x8 -s8 --max-tries=0 --retry-wait=10 -d "$(dirname "$2")" -o "$(basename "$2")" "$1"
      else curl -fL -C - -o "$2" "$1"; fi; }

ena_fastq(){  # $1 = run accession -> prints space-separated ftp urls
  # ENA always returns run_accession as column 1, so fastq_ftp is column 2.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -n +2 | cut -f2 | tr ';' '\n' | sed 's#^#https://#'
}

ONT_FQ=""; R1=""; R2=""
if [ -n "$ONT" ]; then
  urls=$(ena_fastq "$ONT"); u=$(echo "$urls" | head -1)
  ONT_FQ="$OUT/$SAMPLE/${SAMPLE}.ont.fastq.gz"; echo "[sra] ONT $ONT -> $u"; dl "$u" "$ONT_FQ"
  if [ -n "$SUBSAMPLE" ]; then
    command -v seqkit >/dev/null 2>&1 && { seqkit head -n "$SUBSAMPLE" "$ONT_FQ" -o "${ONT_FQ%.gz}.sub.gz" && mv "${ONT_FQ%.gz}.sub.gz" "$ONT_FQ"; }
  fi
fi
if [ -n "$PAIRED" ]; then
  mapfile -t purls < <(ena_fastq "$PAIRED")
  R1="$OUT/$SAMPLE/${SAMPLE}_R1.fastq.gz"; R2="$OUT/$SAMPLE/${SAMPLE}_R2.fastq.gz"
  echo "[sra] Illumina $PAIRED -> ${purls[0]:-} ${purls[1]:-}"
  [ -n "${purls[0]:-}" ] && dl "${purls[0]}" "$R1"
  [ -n "${purls[1]:-}" ] && dl "${purls[1]}" "$R2"
fi

SS="$OUT/samplesheet.$SAMPLE.csv"
echo "sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season" > "$SS"
echo "$SAMPLE,$ONT_FQ,$R1,$R2,$COMP,$FAC,$SEASON" >> "$SS"
echo "[sra] samplesheet -> $SS"
