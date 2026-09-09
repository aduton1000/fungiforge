#!/usr/bin/env bash
# ==============================================================================
# forge audit — where did a *forge HPC install (callforge / captureforge /
# fungiforge) get left off? Read-only. Runs as an ordinary user, no sudo.
#
#   bash hpc_forge_audit.sh                       # audits all three under /hpc/opt
#   bash hpc_forge_audit.sh callforge             # one forge
#   bash hpc_forge_audit.sh callforge:/other/root # explicit install root
#   --out FILE    report path (default $HOME/forge_audit_<host>_<date>.txt)
#
# Per install it reports: repo state (remote, ref, commit, behind origin?, dirty),
# expected layout (launcher, .sif, cli-env, site env, params file, refs), whether
# every path named in the site env / params file exists, PATH integration, the
# effective Nextflow config for the hpc profile, run history (.nextflow/history,
# logs, results) under $HOME and shared project roots, and a DONE/MISSING
# checklist so the next step is obvious.
# ==============================================================================
OUT=""; TARGETS=()
while [ $# -gt 0 ]; do case "$1" in
  --out) OUT="$2"; shift;; -h|--help) sed -n '2,16p' "$0"; exit 0;; *) TARGETS+=("$1");;
esac; shift; done
[ ${#TARGETS[@]} -gt 0 ] || TARGETS=(callforge captureforge fungiforge)
HOST="$(hostname -s 2>/dev/null || hostname)"
[ -n "$OUT" ] || OUT="$HOME/forge_audit_${HOST}_$(date +%Y%m%d_%H%M%S).txt"
TMO=""; command -v timeout >/dev/null 2>&1 && TMO="timeout -k 5"
# Over a non-login ssh (`ssh host 'bash audit.sh'`) /etc/profile.d is NOT sourced, so the
# forge env files and the shared nextflow never reach PATH and every PATH check is a false
# MISSING. Load the site hooks ourselves (they are silent, export-only snippets).
SOURCED_HOOKS=""
for f in /etc/profile.d/*.sh; do
  [ -r "$f" ] && grep -qs 'hpc/opt\|forge' "$f" && { . "$f" >/dev/null 2>&1 || true; SOURCED_HOOKS="$SOURCED_HOOKS $f"; }
done

hdr(){ echo; echo "=================================================================="; echo "## $*"; echo "=================================================================="; }
sub(){ echo; echo "--- $* ---"; }
run(){ local secs="$1" label="$2"; shift 2; echo; echo "### $label"; echo "\$ $*"
  local out rc; if [ -n "$TMO" ]; then out="$($TMO "$secs" bash -c "$*" 2>&1)"; rc=$?; else out="$(bash -c "$*" 2>&1)"; rc=$?; fi
  if [ -n "$out" ]; then printf '%s\n' "$out" | head -300 | sed 's/^/    /'; else echo "    (no output)"; fi
  [ $rc -eq 124 ] && echo "    [TIMED OUT after ${secs}s]"; [ $rc -ne 0 ] && [ $rc -ne 124 ] && echo "    [exit code $rc]"; return 0; }
CHECK=()  # "STATUS|item|detail"
ok(){ CHECK+=("DONE   |$1|$2"); }
miss(){ CHECK+=("MISSING|$1|$2"); }
part(){ CHECK+=("PARTIAL|$1|$2"); }
exists(){ [ -e "$1" ]; }
size(){ du -sh "$1" 2>/dev/null | cut -f1; }

# check every absolute path assigned in a shell env / yaml / nextflow config file
paths_in_file(){ # paths_in_file <file>
  grep -vE '^\s*(#|//)' "$1" 2>/dev/null | grep -oE '(=|:)\s*"?/[A-Za-z0-9_./$@{}-]+' | sed -E 's/^[=:]\s*"?//' | sort -u
}
check_paths(){ # check_paths <file> <label>
  local f="$1" lab="$2" n=0 m=0
  [ -f "$f" ] || return
  sub "paths referenced in $lab ($f)"
  while read -r p; do
    [ -z "$p" ] && continue
    local ep; ep="$(eval echo "$p" 2>/dev/null)"; n=$((n+1))
    if exists "$ep"; then printf '  OK       %-70s %s\n' "$ep" "$(size "$ep")"; else printf '  MISSING  %s\n' "$ep"; m=$((m+1)); fi
  done < <(paths_in_file "$f")
  [ $n -eq 0 ] && echo "  (no absolute paths found)"
  [ $m -gt 0 ] && part "$lab paths" "$m of $n referenced paths missing" || { [ $n -gt 0 ] && ok "$lab paths" "$n paths present"; }
}

audit_forge(){
  local name="$1" root="$2" UP; UP="$(echo "$name" | tr a-z A-Z)"
  CHECK=()
  hdr "$name  —  $root"
  if [ ! -d "$root" ]; then echo "install root does not exist"; miss "install root" "$root"; summary "$name"; return; fi
  run 10 "root" "ls -la '$root'; stat -c '%A %U:%G %n' '$root'"
  # repo: either root itself is the checkout or root/repo
  local repo=""; for c in "$root" "$root/repo"; do [ -f "$c/main.nf" ] && { repo="$c"; break; }; done
  if [ -n "$repo" ]; then
    ok "repo checkout" "$repo"
    run 20 "git state" "cd '$repo' && git remote -v | head -2; git rev-parse --abbrev-ref HEAD; git describe --tags --always 2>/dev/null; git log -1 --format='%h %ci %s'; git status --porcelain | head -10 | sed 's/^/dirty: /'"
    run 40 "behind origin?" "cd '$repo' && git fetch -q origin 2>&1 | head -2; git rev-list --count HEAD..origin/\$(git rev-parse --abbrev-ref HEAD) 2>/dev/null | sed 's/^/commits behind origin: /'; git log -1 --format='origin/main: %h %ci' origin/main 2>/dev/null"
    run 10 "manifest" "grep -E 'version|nextflowVersion' '$repo/nextflow.config' | head -3"
    local dirty; dirty="$(git -C "$repo" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
    [ "$dirty" != "0" ] && part "repo clean" "$dirty modified/untracked files in the shared checkout"
  else miss "repo checkout" "no main.nf under $root or $root/repo"; fi

  sub "layout"
  local f
  for f in "$repo/share/bin/$name-run" "$root/bin/$name-run" "$root/bin/$name" "$repo/.tools/nextflow" "$repo/conf/hpc_slurm.config" "$repo/conf/base.config"; do
    [ -n "$repo" ] || [[ "$f" != "$repo"* ]] || continue
    if [ -e "$f" ]; then printf '  present  %s\n' "$f"; else printf '  absent   %s\n' "$f"; fi
  done
  [ -e "$repo/share/bin/$name-run" ] || [ -e "$root/bin/$name-run" ] && ok "launcher" "$name-run present" || miss "launcher" "share/bin/$name-run"

  sub "container images (.sif) under $root"
  local sifs; sifs="$(find "$root" -maxdepth 3 \( -name '*.sif' -o -name '*.simg' -o -name '*.img' \) -size +1M 2>/dev/null)"
  if [ -n "$sifs" ]; then echo "$sifs" | while read -r s; do printf '  %-8s %s  %s\n' "$(size "$s")" "$(stat -c %y "$s" | cut -d. -f1)" "$s"; done; ok "container image" "$(echo "$sifs" | wc -l | tr -d ' ') image file(s)"
  else echo "  none"; miss "container image" "no .sif under $root (build from env/*.def or docker->SIF)"; fi
  local cache; for cache in "$root/container_cache" "$root/images/cache" "$root/.apptainer" "$root/cache"; do
    [ -d "$cache" ] && run 10 "image cache $cache" "ls -la '$cache' | head -20; du -sh '$cache'"; done
  [ -n "$repo" ] && run 10 "images the pipeline expects (conf/base.config)" "grep -oE \"container *= *[^,}]+\" '$repo/conf/base.config' | sort -u | head -20"
  [ -n "$repo" ] && run 10 "image definitions shipped" "ls '$repo/env/' 2>/dev/null"

  sub "python CLI env"
  local cli; for cli in "$root/cli-env" "$root/env" "$root/venv"; do
    if [ -x "$cli/bin/python" ]; then
      run 20 "cli-env $cli" "'$cli/bin/python' --version; ls '$cli/bin' | grep -iE '^$name' ; '$cli/bin/$name' --version 2>&1 | head -2 || '$cli/bin/$name' version 2>&1 | head -2"
      [ -x "$cli/bin/$name" ] && ok "cli-env" "$cli/bin/$name" || part "cli-env" "$cli exists but no $name entrypoint (pip install the repo into it)"
      break
    fi; done
  [ -x "$root/cli-env/bin/python" ] || [ -x "$root/env/bin/python" ] || [ -x "$root/venv/bin/python" ] || { echo "  none under $root"; miss "cli-env" "no conda env with the $name CLI"; }
  run 10 "conda envs named like $name" "conda env list 2>/dev/null | grep -i '$name'"

  sub "site env / params / config"
  local envf=""; for f in "$root/$name-env.sh" "$repo/share/$name-env.sh" "$HOME/.$name-env.sh" "$repo/$name-env.sh"; do [ -f "$f" ] && { envf="$f"; break; }; done
  if [ -n "$envf" ]; then run 10 "site env $envf" "grep -vE '^\s*#|^\s*$' '$envf'"; ok "site env" "$envf"; check_paths "$envf" "site env"
    grep -q '\[EDIT' "$envf" && part "site env edited" "still contains [EDIT] markers"
  else echo "  no $name-env.sh found (root, repo/share, ~)"; miss "site env" "$name-env.sh not created from share/$name-env.sh.example"; fi
  local params; params="$(find "$root" "$HOME" /hpc/prj /hpc/data 2>/dev/null -maxdepth 3 -name '*.params.yaml' -o -maxdepth 3 -name '*params*.yml' 2>/dev/null | grep -i "$name\|site" | head -5)"
  if [ -n "$params" ]; then for f in $params; do run 10 "params file $f" "grep -vE '^\s*#|^\s*$' '$f' | head -60"; check_paths "$f" "params file"; grep -q '\[EDIT' "$f" && part "params edited" "$f still contains [EDIT] markers"; done; ok "params file" "$params"
  else echo "  no site params yaml found"; [ "$name" != fungiforge ] && miss "params file" "conf/cluster.params.example.yaml not copied/edited"; fi
  for f in "$root/site.config" "$repo/site.config"; do [ -f "$f" ] && { run 10 "site.config" "cat '$f'"; check_paths "$f" "site.config"; ok "site.config" "$f"; }; done

  sub "reference / database roots"
  local refs; refs="$(env | grep -E "^${UP}_(REFS|DB|DATA)=" | cut -d= -f2)"
  for f in $refs "$root/refs" "$root/config" "$root/images" /hpc/refs/$name /hpc/data/$name /hpc/refs; do
    [ -d "$f" ] && run 30 "contents of $f" "ls -la '$f' | head -40; echo; du -sh '$f'/* 2>/dev/null | sort -h | tail -25"; done
  [ -d "$root/refs" ] && run 60 "what kind of refs (fasta / index / vcf / cache)" "cd '$root/refs' && ls | sed -E 's/.*\.(fa|fasta|fna)(\.gz)?$/GENOME &/; s/.*\.(0123|bwt\.2bit\.64|amb|ann|pac|sa)$/BWA-INDEX &/; s/.*\.bt2$/BOWTIE2-INDEX &/; s/.*\.(fai|dict)$/GENOME-AUX &/; s/.*\.vcf(\.gz)?$/VCF &/; s/.*\.(bw|bigwig)$/BIGWIG &/; s/.*\.bed(\.gz)?$/BED &/; s/.*\.gtf(\.gz)?$/GTF &/' | sort | head -60; find . -maxdepth 2 -type d | head -20"
  [ -n "$refs" ] && { for f in $refs; do [ -d "$f" ] && ok "refs root" "$f ($(size "$f"))" || miss "refs root" "$f (named in env) does not exist"; done; }

  sub "PATH integration"
  run 10 "commands on PATH" "command -v $name $name-run; type $name-run 2>/dev/null | head -2"
  run 10 "who puts $root on PATH" "grep -ls '$root' /etc/profile.d/* /etc/profile /etc/bash.bashrc /etc/environment \$HOME/.bashrc \$HOME/.profile \$HOME/.bash_profile 2>/dev/null; grep -h '$root' /etc/profile.d/* \$HOME/.bashrc 2>/dev/null | head -5"
  command -v "$name-run" >/dev/null 2>&1 && ok "PATH" "$name-run resolves to $(command -v "$name-run")" || miss "PATH" "$name-run not on PATH (profile.d snippet / env file not sourced)"
  run 10 "$UP_* env vars" "env | grep -E '^${UP}_' | sort"

  if [ -n "$repo" ] && command -v nextflow >/dev/null 2>&1; then
    sub "effective nextflow config (hpc profile)"
    local eng=apptainer; command -v apptainer >/dev/null 2>&1 || eng=singularity
    run 90 "nextflow config -profile hpc_slurm,$eng" "cd '$repo' && nextflow config -profile hpc_slurm,$eng 2>&1 | grep -E \"^(process|executor|params\.(slurm|max_|container|bind|scratch|data_dir|.*image)|apptainer|singularity)\" | head -60"
  fi

  sub "run history (where did we leave off?)"
  run 60 "nextflow history files mentioning $name" "for h in \$(find \$HOME /hpc/prj /hpc/data -maxdepth 4 -path '*/.nextflow/history' 2>/dev/null); do if grep -qi '$name' \"\$h\" 2>/dev/null; then echo \"== \$h\"; tail -5 \"\$h\" | cut -c1-200; fi; done"
  run 60 "recent nextflow logs mentioning $name (last status line)" "for l in \$(find \$HOME /hpc/prj -maxdepth 3 -name '.nextflow.log*' -newermt '-120 days' 2>/dev/null | head -30); do grep -qi '$name' \"\$l\" 2>/dev/null || continue; echo \"== \$l (\$(stat -c %y \"\$l\" | cut -d. -f1))\"; grep -E 'Launching|Execution complete|Execution cancelled|Session aborted|ERROR ~|Goodbye|Pipeline completed' \"\$l\" | tail -3 | cut -c1-220; done"
  run 30 "last run's launch command + params (from history)" "for h in \$(find \$HOME /hpc/prj -maxdepth 4 -path '*/.nextflow/history' 2>/dev/null); do grep -i '$name' \"\$h\" 2>/dev/null | tail -1 | cut -f7- | cut -c1-600; done; ls -d \$HOME/*/results* 2>/dev/null | head; for r in \$(ls -d \$HOME/cf_sim*/results* \$HOME/cf_smoke/results* 2>/dev/null | head -4); do echo \"== \$r\"; ls \"\$r\" | head -12; done"
  run 30 "results / work dirs mentioning $name" "find \$HOME /hpc/prj -maxdepth 3 -type d \( -iname 'results*' -o -iname '*${name}*' -o -name 'work' \) -newermt '-365 days' 2>/dev/null | head -20"
  run 20 "handover / notes in the install" "ls '$root'/*.md '$repo'/HANDOVER.md '$repo'/docs/*.md 2>/dev/null; head -30 '$repo/HANDOVER.md' 2>/dev/null"
  summary "$name"
}

summary(){
  local name="$1"
  sub "CHECKLIST — $name"
  local c; for c in "${CHECK[@]}"; do IFS='|' read -r st item det <<<"$c"; printf '  %-8s %-22s %s\n' "$st" "$item" "$det"; done
  local nm; nm="$(printf '%s\n' "${CHECK[@]}" | grep -c '^MISSING')"; local np; np="$(printf '%s\n' "${CHECK[@]}" | grep -c '^PARTIAL')"
  echo; echo "  => $name: $nm missing, $np partial"
}

{
  echo "forge audit — $(date) — $(id -un)@$(hostname -f 2>/dev/null || hostname)"
  echo "targets: ${TARGETS[*]}"
  hdr "shared context"
  echo "profile.d hooks sourced by this audit:${SOURCED_HOOKS:- none}"
  run 10 "PATH" "echo \"\$PATH\" | tr ':' '\n' | awk '!seen[\$0]++' | grep -E 'hpc|opt|forge'"
  run 10 "/hpc/opt contents" "ls -la /hpc/opt 2>/dev/null; ls -la /hpc/opt/*/ 2>/dev/null | head -60"
  run 10 "profile.d hooks" "ls -la /etc/profile.d/ 2>/dev/null; grep -l 'hpc/opt' /etc/profile.d/* 2>/dev/null | xargs -r cat"
  run 10 "toolchain" "nextflow -v 2>&1 | head -1; apptainer --version 2>&1; singularity --version 2>&1 | head -1; docker --version 2>&1; java -version 2>&1 | head -1; conda --version 2>&1; mamba --version 2>&1 | head -1"
  run 20 "docker images on this host (forge-related)" "docker images --format '{{.Repository}}:{{.Tag}}  {{.Size}}  {{.CreatedAt}}' 2>/dev/null | grep -iE 'forge|antismash|funannotate|flye|medaka|busco|tetools|vep|deepvariant|gatk|hap' | head -30"
  run 10 "apptainer caches" "du -sh \$HOME/.apptainer \$HOME/.singularity 2>/dev/null; ls \$HOME/.apptainer/cache 2>/dev/null"
  run 10 "/hpc/refs /hpc/data top level" "ls -la /hpc/refs /hpc/data 2>/dev/null | head -40"
  run 10 "ownership of shared roots (who can install?)" "stat -c '%A %U:%G %n' /hpc /hpc/opt /hpc/opt/* /hpc/data /hpc/data/* /hpc/prj /hpc/refs 2>/dev/null; id"
  run 10 "docker daemon access (image builds)" "docker ps >/dev/null 2>&1 && echo DOCKER_OK || echo 'no docker access'; docker system df 2>/dev/null | head -5"
  run 10 "existing kraken2 DBs (reusable for decontam)" "ls -la /hpc/data/kraken2 /data/kraken2 2>/dev/null | head -20; env | grep -i kraken"
  run 10 "how /hpc/opt/*/bin lands on PATH" "grep -rn 'hpc/opt\|hpc/bin' /etc/profile.d/ /etc/profile /etc/bash.bashrc /etc/environment /hpc/etc 2>/dev/null | head -20; ls -la /hpc/bin /hpc/sbin 2>/dev/null | head -20"
  for t in "${TARGETS[@]}"; do
    name="${t%%:*}"; root="${t#*:}"; [ "$root" = "$t" ] && root="/hpc/opt/$name"
    audit_forge "$name" "$root"
  done
  hdr "END"
  echo "Bring this file back: $OUT"
} 2>&1 | tee "$OUT"
echo; echo "Report written: $OUT ($(wc -l < "$OUT") lines)"
