#!/usr/bin/env bash
# Run the three-mode read fixture through the REAL containers with Docker (not -stub-run).
#
# What it checks that the stub DAG cannot: the images resolve and start, bin/ is on PATH inside
# every container (`source ff_status.sh`, `check_databases.py`), the status contract records
# versions and exit codes in-container, and DB_CHECK behaves without databases
# (--allow_missing_db). The fixture reads are tiny synthetic files, so the assemblers are
# EXPECTED to fail on them: the run stops at ASSEMBLE / SR_ASSEMBLE with `failed` recorded in
# the assemble JSONs. Success here means: DB_CHECK + every READ_QC task completed, and the
# assembler failures were recorded by the contract rather than hidden.
#
# Needs: docker, nextflow, the images from conf/base.config pulled locally (or pullable).
# This is a local/dev-deployment check; it is not part of the hosted CI (image set ~30 GB).
#
#   bash test/run_docker_fixture.sh [work_dir]
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$REPO/work_docker_fixture}"
OUT="$REPO/results_docker_fixture"
rm -rf "$OUT"
mkdir -p "$REPO/test/mini_db"
cd "$REPO" || exit 1
nextflow run main.nf -profile local,test,docker --outdir "$OUT" --allow_missing_db -w "$WORK" > "$OUT.log" 2>&1
rc=$?
echo "nextflow exit: $rc (non-zero is expected: the fixture reads cannot be assembled)"

ok=1
check() { if eval "$2"; then echo "  ok    $1"; else echo "  FAIL  $1"; ok=0; fi; }
check "DB_CHECK ran"                 "grep -q 'DB_CHECK' '$OUT.log'"
check "READ_QC ran for 3 isolates"   "[ \"\$(grep -c 'READ_QC (TESTFUN0' '$OUT.log')\" -eq 3 ]"
check "readqc JSONs carry versions"  "[ \"\$(find '$WORK' -name '*.readqc.json' -exec grep -l '\"versions\"' {} + | wc -l)\" -ge 3 ]"
check "assembler failure recorded"   "grep -q '\[ff:assemble\].*failed' '$OUT.log'"
check "no silent success"            "! grep -q 'SUCCESS' '$OUT.log'"
if [ "$ok" = 1 ]; then echo "docker fixture: PASS"; exit 0; else echo "docker fixture: FAIL (see $OUT.log)"; exit 1; fi
