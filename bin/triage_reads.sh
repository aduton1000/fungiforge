#!/usr/bin/env bash
# ==============================================================================
# fungiforge — read-level triage: what is actually in each sample, BEFORE assembly
# ------------------------------------------------------------------------------
# Classifies a subsample of reads per isolate with Kraken2 (the PlusPF database
# fungiforge ships in $FUNGIFORGE_DB/kraken2) and reports the domain composition,
# top species and a verdict: fungal | likely_fungal | non_fungal | human | mixed.
# Seconds per sample, so run it on a whole plate before spending hours assembling
# isolates that turn out to be bacterial, human-contaminated, or mixed cultures.
#
# NB the PlusPF-8 database holds only a small set of fungal genomes: a genuine fungal
# isolate typically reads as mostly *unclassified* with a fungal minority ("likely_fungal"),
# whereas a bacterial isolate classifies overwhelmingly as Bacteria.
#
# Usage:
#   triage_reads.sh --samplesheet samples.csv [--out triage] [--n 200000] [--threads 8]
#   triage_reads.sh --r1 'reads/*_R1_001.fastq.gz' --r2 'reads/*_R2_001.fastq.gz' [...]
#   triage_reads.sh --ont 'reads/*.fastq.gz' [...]
# Needs kraken2 + python3 on PATH (both are in the fungiforge image; on a cluster:
#   apptainer exec -B /hpc <fungiforge.sif> bash bin/triage_reads.sh ...), and the DB:
#   --db DIR  (default $FUNGIFORGE_DB/kraken2).
# Outputs in --out: <sample>.k2.report, triage.tsv, and — when --samplesheet is given —
#   samples.fungal.csv (fungal + likely_fungal rows, plus mixed rows with more fungal than
#   bacterial reads — Stage 04 decontaminates those) and samples.excluded.csv.
# ==============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHEET=""; R1=(); R2=(); ONT=(); OUT=triage; N=200000; THREADS=8; DB="${FUNGIFORGE_DB:+$FUNGIFORGE_DB/kraken2}"
while [ $# -gt 0 ]; do case "$1" in
  --samplesheet) SHEET="$2"; shift;; --out) OUT="$2"; shift;; --n) N="$2"; shift;; --threads) THREADS="$2"; shift;;
  --db) DB="$2"; shift;;
  --r1) shift; while [ $# -gt 0 ] && [[ "$1" != --* ]]; do R1+=("$1"); shift; done; continue;;
  --r2) shift; while [ $# -gt 0 ] && [[ "$1" != --* ]]; do R2+=("$1"); shift; done; continue;;
  --ont) shift; while [ $# -gt 0 ] && [[ "$1" != --* ]]; do ONT+=("$1"); shift; done; continue;;
  -h|--help) sed -n '2,26p' "$0"; exit 0;; *) echo "unknown option $1" >&2; exit 2;;
esac; shift; done
[ -n "$DB" ] && [ -f "$DB/hash.k2d" ] || { echo "ERROR: Kraken2 DB not found (--db or \$FUNGIFORGE_DB/kraken2): '$DB'" >&2; exit 2; }
command -v kraken2 >/dev/null || { echo "ERROR: kraken2 not on PATH (run inside the fungiforge image)" >&2; exit 2; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 2; }
mkdir -p "$OUT"; TAX="$HERE/kraken_taxonomy.py"

# ---- collect samples: id, mode, files -------------------------------------
declare -a SID SMODE SF1 SF2
add(){ SID+=("$1"); SMODE+=("$2"); SF1+=("$3"); SF2+=("${4:-}"); }
if [ -n "$SHEET" ]; then
  base="$(cd "$(dirname "$SHEET")" && pwd)"
  resolve(){ [ -z "$1" ] && { echo ""; return; }; case "$1" in /*) echo "$1";; *) [ -e "$1" ] && echo "$PWD/$1" || echo "$base/$1";; esac; }
  tail -n +2 "$SHEET" | while IFS=, read -r s ont r1 r2 rest; do
    [ -z "$s" ] && continue
    if [ -n "$r1" ] && [ -n "$r2" ]; then echo "$s	illumina	$(resolve "$r1")	$(resolve "$r2")"
    elif [ -n "$ont" ]; then echo "$s	ont	$(resolve "$ont")	"
    fi
  done > "$OUT/.samples.tsv"
  while IFS=$'\t' read -r s m f1 f2; do add "$s" "$m" "$f1" "$f2"; done < "$OUT/.samples.tsv"
else
  ilmn_sid(){ local b; b="$(basename "$1")"; b="${b%.fastq.gz}"; b="${b%.fq.gz}"; b="${b%.fastq}"; b="${b%.fq}"
              b="$(echo "$b" | sed -E 's/_S[0-9]+(_L[0-9]{3})?_R[12]_[0-9]{3}$//; s/[._]R?[12]$//')"; echo "$b"; }
  declare -A M2
  for pat in "${R2[@]}"; do for f in $pat; do M2["$(ilmn_sid "$f")"]="$f"; done; done
  for pat in "${R1[@]}"; do for f in $pat; do s="$(ilmn_sid "$f")"; add "$s" illumina "$f" "${M2[$s]:-}"; done; done
  for pat in "${ONT[@]}"; do for f in $pat; do b="$(basename "$f")"; add "${b%%.*}" ont "$f" ""; done; done
fi
[ ${#SID[@]} -gt 0 ] || { echo "no samples found" >&2; exit 2; }
echo "[triage] ${#SID[@]} samples · $N reads(pairs) each · DB $DB · threads $THREADS" >&2

# ---- classify a subsample per sample ----------------------------------------
printf 'sample\tmode\treads_classified_from\tpct_unclassified\tpct_bacteria\tpct_archaea\tpct_viruses\tpct_fungi\tpct_other_eukaryota\tpct_human\ttop_species\ttop_species_pct\tverdict\n' > "$OUT/triage.tsv"
sub(){ gzip -dc "$1" 2>/dev/null | head -n $((N * 4)); }   # first N reads; gz or plain
for i in "${!SID[@]}"; do
  s="${SID[$i]}"; m="${SMODE[$i]}"; f1="${SF1[$i]}"; f2="${SF2[$i]}"; rep="$OUT/$s.k2.report"
  [ -s "$f1" ] || { echo "[triage] $s: missing $f1 — skipped" >&2; continue; }
  if [ "$m" = illumina ] && [ -s "$f2" ]; then
    sub "$f1" > "$OUT/.$s.1.fq"; sub "$f2" > "$OUT/.$s.2.fq"
    kraken2 --db "$DB" --threads "$THREADS" --paired --report "$rep" --output /dev/null "$OUT/.$s.1.fq" "$OUT/.$s.2.fq" >/dev/null 2>"$OUT/.$s.k2.log"
  else
    sub "$f1" > "$OUT/.$s.1.fq"
    kraken2 --db "$DB" --threads "$THREADS" --report "$rep" --output /dev/null "$OUT/.$s.1.fq" >/dev/null 2>"$OUT/.$s.k2.log"
  fi
  rm -f "$OUT/.$s".[12].fq
  if [ -s "$rep" ]; then
    line="$(python3 "$TAX" summarize "$rep")"
    printf '%s\t%s\t%s\n' "$s" "$m" "$line" >> "$OUT/triage.tsv"
    echo "[triage] $s: $(echo "$line" | awk -F'\t' '{printf "%s (%s%% unclassified, %s%% bacteria, %s%% fungi; top %s %s%%)", $11, $2, $3, $6, $9, $10}')" >&2
  else
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tkraken2_failed\n' "$s" "$m" >> "$OUT/triage.tsv"
    echo "[triage] $s: kraken2 failed — see $OUT/.$s.k2.log" >&2
  fi
done

# ---- summary + filtered samplesheets ------------------------------------------
echo >&2; echo "[triage] verdict counts:" >&2; tail -n +2 "$OUT/triage.tsv" | cut -f13 | sort | uniq -c | sed 's/^/    /' >&2
if [ -n "$SHEET" ]; then
  head -1 "$SHEET" > "$OUT/samples.fungal.csv"; head -1 "$SHEET" > "$OUT/samples.excluded.csv"
  # fungal + likely_fungal, plus "mixed" samples where fungal reads outnumber bacterial ones —
  # a contaminated fungal culture is exactly what Stage 04 decontamination is for.
  awk -F'\t' 'NR>1 && ($13=="fungal" || $13=="likely_fungal" || ($13=="mixed" && $8+0 > $5+0)){print $1}' "$OUT/triage.tsv" | sort -u > "$OUT/.keep"
  tail -n +2 "$SHEET" | while IFS= read -r row; do s="${row%%,*}"; if grep -qx "$s" "$OUT/.keep"; then echo "$row" >> "$OUT/samples.fungal.csv"; else echo "$row" >> "$OUT/samples.excluded.csv"; fi; done
  echo "[triage] $(($(wc -l < "$OUT/samples.fungal.csv")-1)) fungal/likely_fungal -> $OUT/samples.fungal.csv ; $(($(wc -l < "$OUT/samples.excluded.csv")-1)) excluded -> $OUT/samples.excluded.csv" >&2
  rm -f "$OUT/.keep" "$OUT/.samples.tsv"
fi
echo "[triage] table: $OUT/triage.tsv" >&2
