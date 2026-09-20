#!/usr/bin/env bash
# Stream a possibly-compressed file to stdout, choosing the decompressor by MAGIC BYTES rather
# than by file name. A release fetched from a mirror can carry one compression under another
# extension — RVDB-prot was staged as rvdb.fasta.gz from a .xz URL — and `gzip -dc` then fails
# with "not in gzip format", which reads like a corrupt download rather than a naming mismatch.
#
#   smart_cat.sh <file>          # gzip | xz | bzip2 | zstd | plain, decided by content
set -uo pipefail
f="${1:?usage: smart_cat.sh <file>}"
[ -r "$f" ] || { echo "smart_cat: cannot read '$f'" >&2; exit 1; }
magic="$(head -c 6 "$f" 2>/dev/null | od -An -tx1 | tr -d ' \n')"
case "$magic" in
  1f8b*)       exec gzip  -dc "$f" ;;
  fd377a585a*) exec xz    -dc "$f" ;;
  425a68*)     exec bzip2 -dc "$f" ;;
  28b52ffd*)   exec zstd  -dc "$f" ;;
  *)           exec cat       "$f" ;;
esac
