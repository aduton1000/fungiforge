#!/usr/bin/env bash
# Run `mlst` against the fungal PubMLST schemes staged by bin/fetch_references.sh mlst (W2.3).
#   mlst_fungal.sh <data_dir>/mlst/pubmlst <assembly.fa> <threads> <out.tsv>
# mlst ships bacterial schemes only, so its BLAST database is built here, in the task
# directory, from the staged <scheme>/<locus>.tfa files (headers become >scheme.locus_allele,
# the form mlst-make_blast_db produces), then mlst is pointed at both with --datadir/--blastdb.
set -euo pipefail
PUBMLST="$1"; ASM="$2"; THREADS="${3:-2}"; OUT="${4:-mlst.tsv}"
mkdir -p mlstdb
: > mlstdb/mlst.fa
for d in "$PUBMLST"/*/; do
  s=$(basename "$d")
  for f in "$d"/*.tfa; do
    [ -e "$f" ] || continue
    sed "s/^>/>$s./" "$f" >> mlstdb/mlst.fa
  done
done
[ -s mlstdb/mlst.fa ] || { echo "no allele files under $PUBMLST" >&2; exit 2; }
makeblastdb -hash_index -in mlstdb/mlst.fa -dbtype nucl -parse_seqids > mlstdb/makeblastdb.log 2>&1
mlst --datadir "$PUBMLST" --blastdb mlstdb/mlst.fa --threads "$THREADS" "$ASM" > "$OUT"
