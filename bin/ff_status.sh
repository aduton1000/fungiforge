#!/usr/bin/env bash
# fungiforge — stage status contract (bash side). Source at the top of every task script:
#
#   source ff_status.sh          # bin/ is on PATH inside every task
#   ff_init  <stage> <sample> <stage.json> [--best-effort]
#   ff_run   <label> [--optional] -- <command ...>      # run a tool, record its exit code (in $FF_RC)
#   ff_skip  <label> <reason>                            # record a tool deliberately not run
#   ff_version <label> -- <command printing a version>   # record the tool's version string
#   ff_finalize                                           # write status into <stage.json>; exit code
#
# Semantics
#   * A REQUIRED tool (default) that exits non-zero marks the stage `failed`; ff_run finalizes
#     immediately and exits — non-zero for a normal stage (the task fails, Nextflow stops the
#     run after pending tasks), zero for a --best-effort stage (status recorded, run continues).
#   * An --optional tool that exits non-zero marks the stage `partial` and execution continues
#     (the module must have a documented fallback for that step). ff_run itself always returns
#     0 in that case — task scripts run under `bash -ue`, so a non-zero return from a bare
#     ff_run statement would abort the task. Test $FF_RC to branch on the tool's outcome.
#   * ff_skip records a tool that was not run and why (missing database, not applicable);
#     it does not change the status.
#   * Every tool's exit code and wall time land in the stage JSON under "tools", so nothing
#     that ran is invisible in the results. `|| true` is never needed and must not be used.
#   * ff_version runs a cheap version command and records the first line that looks like a
#     version (else the first line of output; "not found"/"unknown" when nothing usable).
#     Versions land under "versions" and are copied onto the matching "tools" entry, so the
#     provenance stage can build one version table for the whole run.
#
# Pipelines: wrap in bash with pipefail so the exit code is the pipeline's, e.g.
#   ff_run chopper -- bash -o pipefail -c "chopper -q 10 -i in.fq | gzip > out.fq.gz"

ff_init() {
  FF_STAGE="$1"; FF_SAMPLE="$2"; FF_JSON="$3"; FF_BEST_EFFORT=""
  [ "${4:-}" = "--best-effort" ] && FF_BEST_EFFORT="--best-effort"
  FF_TOOLS=".ff_tools.tsv"; FF_SKIPS=".ff_skips.tsv"; FF_VERSIONS=".ff_versions.tsv"; FF_RC=0
  : > "$FF_TOOLS"; : > "$FF_SKIPS"; : > "$FF_VERSIONS"
  FF_HELPER="${FF_HELPER:-$(dirname "${BASH_SOURCE[0]}")/ff_status.py}"
}

ff_version() {
  local label="$1"; shift
  [ "${1:-}" = "--" ] && shift
  local out v
  if ! command -v "$1" >/dev/null 2>&1; then
    v="not found"
  else
    # never let a broken --version abort the task: capture, cap, pick the version-looking line
    out=$({ "$@" 2>&1 || true; } | head -c 4000 | tr -d '\r')
    v=$(printf '%s\n' "$out" | grep -m1 -E '[0-9]+\.[0-9]+' || true)
    [ -z "$v" ] && v=$(printf '%s\n' "$out" | grep -m1 . || true)
    [ -z "$v" ] && v="unknown"
    v=$(printf '%s' "$v" | tr -s '[:space:]' ' ' | sed -E 's/^ +| +$//g' | cut -c1-200)
  fi
  printf '%s\t%s\n' "$label" "$v" >> "$FF_VERSIONS"
  return 0
}

ff_run() {
  local label="$1"; shift
  local optional=0
  while [ $# -gt 0 ]; do
    case "$1" in
      --optional) optional=1; shift ;;
      --) shift; break ;;
      *) break ;;
    esac
  done
  local t0=$SECONDS rc=0
  # `cmd && rc=0 || rc=$?` keeps errexit from aborting inside the function
  "$@" && rc=0 || rc=$?
  # shellcheck disable=SC2034  # FF_RC is read by the calling task script
  FF_RC=$rc
  printf '%s\t%s\t%s\t%s\n' "$label" "$rc" "$optional" "$((SECONDS - t0))" >> "$FF_TOOLS"
  if [ "$rc" -ne 0 ]; then
    if [ "$optional" -eq 1 ]; then
      echo "[ff:$FF_STAGE] $label exited $rc (optional step; continuing with the documented fallback)" >&2
    else
      echo "[ff:$FF_STAGE] $label exited $rc (required step) — stage failed" >&2
      ff_finalize
      exit $?
    fi
  fi
  return 0
}

ff_skip() {
  printf '%s\t%s\n' "$1" "$2" >> "$FF_SKIPS"
  echo "[ff:$FF_STAGE] $1 skipped: $2" >&2
}

ff_finalize() {
  python3 "$FF_HELPER" finalize --stage "$FF_STAGE" --sample "$FF_SAMPLE" --json "$FF_JSON" \
      --tools "$FF_TOOLS" --skips "$FF_SKIPS" --versions "${FF_VERSIONS:-}" $FF_BEST_EFFORT
}
