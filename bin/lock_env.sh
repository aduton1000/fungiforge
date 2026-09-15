#!/usr/bin/env bash
# Regenerate env/base.linux-64.lock — the explicit (URL + md5) lock of the fungiforge base
# environment — from a built image, and check that env/base.yml's exact pins solve.
#
#   bin/lock_env.sh [image]        default image: aduton1000/fungiforge:<version in nextflow.config>
#
# Workflow when upgrading a tool: edit the `==` pin in env/base.yml -> build the image from
# base.yml once (docker build --build-arg ENV_SPEC=base.yml, or temporarily point the
# Dockerfile at base.yml) -> validate the run -> bin/lock_env.sh <image> -> commit the lock.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VER=$(grep -oE "version *= *'[^']+'" "$REPO/nextflow.config" | head -1 | sed -E "s/.*'([^']+)'/\1/")
IMG="${1:-aduton1000/fungiforge:$VER}"
OUT="$REPO/env/base.linux-64.lock"
command -v docker >/dev/null || { echo "docker not found" >&2; exit 1; }
echo "[lock_env] exporting explicit lock from $IMG"
{
  echo "# fungiforge base environment — explicit lock exported from the validated image"
  echo "# $IMG (micromamba env export --explicit). Regenerate: bin/lock_env.sh"
  docker run --rm --platform linux/amd64 "$IMG" bash -c 'micromamba env export -n base --explicit' \
    | grep -v '^# This file may\|^# \$ conda'
} > "$OUT.tmp"
grep -q '^@EXPLICIT' "$OUT.tmp" || { echo "[lock_env] export failed" >&2; rm -f "$OUT.tmp"; exit 1; }
mv "$OUT.tmp" "$OUT"
echo "[lock_env] $(grep -c '^https' "$OUT") packages -> $OUT"
echo "[lock_env] checking that env/base.yml's pins agree with the lock"
python3 - "$REPO/env/base.yml" "$OUT" <<'PY'
import re, sys
yml, lock = sys.argv[1], sys.argv[2]
pins = dict(re.findall(r'^\s*- ([A-Za-z0-9_-]+)==([^\s#]+)', open(yml).read(), re.M))
have = {}
for line in open(lock):
    if line.startswith("https"):
        fn = line.rsplit("/", 1)[1].split("#")[0]
        m = re.match(r'(.+?)-(\d[^-]*)-[^-]+\.(conda|tar\.bz2)$', fn)
        if m: have[m.group(1).lower()] = m.group(2)
bad = {k: (v, have.get(k)) for k, v in pins.items() if have.get(k.lower()) != v}
if bad:
    print("[lock_env] MISMATCH between base.yml pins and the image:", bad); sys.exit(1)
print(f"[lock_env] all {len(pins)} pins match the lock")
PY
