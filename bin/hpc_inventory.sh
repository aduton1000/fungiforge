#!/usr/bin/env bash
# ==============================================================================
# fungiforge — HPC inventory / pre-installation site survey
# ------------------------------------------------------------------------------
# Read-only survey of a cluster to plan a site-wide fungiforge deployment
# (SLURM + Apptainer/Singularity + Nextflow + Lmod). Runs as an ordinary user,
# needs no sudo, changes nothing on the system. Writes ONE text report to $HOME
# that you bring back for review.
#
# Usage (on an HPC login node):
#   bash hpc_inventory.sh                      # login-node survey only
#   bash hpc_inventory.sh --submit-test        # + submit a 10-min SLURM job that
#                                              #   surveys a compute node too
#   bash hpc_inventory.sh --pull-test          # + pull a 3 MB test image with
#                                              #   apptainer (proves Docker Hub access)
#   bash hpc_inventory.sh --submit-test --pull-test --partition short --account myacct
#
# Options:
#   --submit-test          submit a short SLURM job running the compute-node probe
#   --pull-test            apptainer pull docker://alpine into a temp cache dir
#   --partition <name>     partition for the test job (default: cluster default)
#   --account <name>       account for the test job
#   --wait <seconds>       how long to wait for the test job (default 900)
#   --out <file>           report path (default $HOME/fungiforge_hpc_inventory_<host>_<date>.txt)
#   --db-root <dir>        extra directory to scan for existing reference databases
#
# What it checks: OS/arch/kernel · scheduler (partitions, limits, accounts, QoS)
# · module system + available modules · nextflow/java/apptainer/singularity/
# docker/podman/conda/aria2c on PATH · user-namespace + setuid container support
# · storage (mounts, free space, quotas, writable shared install roots, spaces in
# paths) · outbound network to every host fungiforge downloads from · existing
# bioinformatics databases (BUSCO, funannotate, antiSMASH, eggNOG, Kraken2,
# UNITE, InterProScan, GeneMark key) · existing fungiforge installs · login-node
# resource limits (ulimit, cgroup) that affect a Nextflow head process.
# ==============================================================================

INVENTORY_VERSION="1.0"
MODE="login"
SUBMIT_TEST=0; PULL_TEST=0; PARTITION=""; ACCOUNT=""; WAIT_SECS=900; OUT=""; EXTRA_DB_ROOT=""

while [ $# -gt 0 ]; do
  case "$1" in
    --compute-probe) MODE="compute" ;;
    --submit-test)   SUBMIT_TEST=1 ;;
    --pull-test)     PULL_TEST=1 ;;
    --partition)     PARTITION="$2"; shift ;;
    --account)       ACCOUNT="$2"; shift ;;
    --wait)          WAIT_SECS="$2"; shift ;;
    --out)           OUT="$2"; shift ;;
    --db-root)       EXTRA_DB_ROOT="$2"; shift ;;
    -h|--help)       sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
HOST_SHORT="$(hostname -s 2>/dev/null || hostname)"
STAMP="$(date +%Y%m%d_%H%M%S)"
[ -n "$OUT" ] || OUT="$HOME/fungiforge_hpc_inventory_${HOST_SHORT}_${STAMP}.txt"
TMO=""; command -v timeout >/dev/null 2>&1 && TMO="timeout -k 5"

# ---- helpers -----------------------------------------------------------------
hdr(){ echo; echo "=================================================================="; echo "## $*"; echo "=================================================================="; }
sub(){ echo; echo "--- $* ---"; }
# run <seconds> <label> <shell command string>   (runs in a child bash, timed, indented)
run(){
  local secs="$1" label="$2"; shift 2; local cmd="$*"
  echo; echo "### $label"; echo "\$ $cmd"
  local out rc
  if [ -n "$TMO" ]; then out="$($TMO "$secs" bash -c "$cmd" 2>&1)"; rc=$?
  else out="$(bash -c "$cmd" 2>&1)"; rc=$?; fi
  if [ -n "$out" ]; then printf '%s\n' "$out" | head -400 | sed 's/^/    /'; else echo "    (no output)"; fi
  [ $rc -eq 124 ] && echo "    [TIMED OUT after ${secs}s]"
  [ $rc -ne 0 ] && [ $rc -ne 124 ] && echo "    [exit code $rc]"
  return 0
}
# has <cmd>
has(){ command -v "$1" >/dev/null 2>&1; }
# tool <cmd> <version-args>: one-line presence + version + path
tool(){
  local c="$1"; shift
  if has "$c"; then
    local v; v="$($TMO 20 "$c" "$@" 2>&1 | head -3 | tr '\n' ' ' | cut -c1-160)"
    printf '  %-14s FOUND  %s\n' "$c" "$(command -v "$c")"
    printf '  %-14s        %s\n' "" "$v"
  else
    printf '  %-14s missing\n' "$c"
  fi
}
# init the module system in THIS shell so `module` works in eval'd commands
init_modules(){
  if ! type module >/dev/null 2>&1; then
    for f in "$LMOD_PKG/init/bash" "$MODULESHOME/init/bash" /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
             /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash /opt/apps/lmod/lmod/init/bash \
             /cm/local/apps/environment-modules/current/init/bash; do
      [ -n "$f" ] && [ -r "$f" ] && { source "$f" >/dev/null 2>&1; type module >/dev/null 2>&1 && break; }
    done
  fi
  type module >/dev/null 2>&1
}
# runm <label> <module command...>  (runs in current shell; module output is on stderr)
runm(){
  local label="$1"; shift
  echo; echo "### $label"; echo "\$ $*"
  local out; out="$(eval "$*" 2>&1)"
  if [ -n "$out" ]; then printf '%s\n' "$out" | head -300 | sed 's/^/    /'; else echo "    (no output)"; fi
}
# url_check <url>: http code + timing, honours proxy env
url_check(){
  local url="$1" code
  code="$($TMO 25 curl -sS -L -o /dev/null -m 20 -w '%{http_code} %{time_total}s %{remote_ip}' "$url" 2>&1 | tail -1)"
  printf '  %-58s %s\n' "$url" "$code"
}

# ==============================================================================
# COMPUTE-NODE PROBE (also runs inside the SLURM test job)
# ==============================================================================
compute_probe(){
  hdr "COMPUTE-NODE PROBE  (host: $(hostname)  job: ${SLURM_JOB_ID:-none}  $(date))"
  run 10 "node identity"      "hostname -f; uname -srm; cat /etc/os-release | head -4"
  run 10 "cpu / memory"       "nproc; free -g; lscpu | grep -E 'Model name|Socket|Thread|^CPU\(s\)|Flags' | cut -c1-200"
  run 10 "SLURM job env"      "env | grep -E '^SLURM_(JOB_ID|JOB_PARTITION|JOB_ACCOUNT|CPUS_ON_NODE|MEM_PER_NODE|JOB_NODELIST|TMPDIR)' ; echo TMPDIR=\${TMPDIR:-unset}"
  run 10 "cgroup limits"      "cat /sys/fs/cgroup/memory.max 2>/dev/null; cat /sys/fs/cgroup/cpu.max 2>/dev/null; cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null"
  run 10 "gpu"                "nvidia-smi -L 2>&1 | head -5"
  run 20 "mounts visible"     "df -hP \$HOME /tmp /scratch /work /project /projects /data /shared /apps /opt /cvmfs 2>/dev/null"
  run 20 "filesystem types"   "mount | awk '{print \$5, \$3}' | grep -Ei 'lustre|gpfs|nfs|beegfs|ceph|xfs|ext4|tmpfs|overlay' | sort | uniq -c | sort -rn | head -30"
  run 10 "local tmp size"     "df -hP /tmp \${TMPDIR:-/tmp} 2>/dev/null; ls -ld /tmp \${TMPDIR:-/tmp}"
  run 10 "user namespaces"    "echo max_user_namespaces=\$(cat /proc/sys/user/max_user_namespaces 2>/dev/null); echo unprivileged_userns_clone=\$(cat /proc/sys/kernel/unprivileged_userns_clone 2>/dev/null); grep -E \"^\$(id -un):\" /etc/subuid /etc/subgid 2>/dev/null"
  sub "container runtimes on compute node"
  tool apptainer --version; tool singularity --version; tool docker --version; tool podman --version
  tool nextflow -v; tool java -version; tool aria2c --version; tool curl --version
  if init_modules; then
    runm "modules on compute node" "module --terse avail 2>&1 | grep -iE 'nextflow|java|jdk|openjdk|apptainer|singularity|conda|mamba|python' | head -40"
  else
    echo "    (no module system found on compute node)"
  fi
  sub "outbound network from compute node"
  for u in https://github.com https://registry-1.docker.io/v2/ https://quay.io https://conda.anaconda.org https://get.nextflow.io \
           https://busco-data.ezlab.org https://unite.ut.ee https://genome-idx.s3.amazonaws.com https://ftp.ncbi.nlm.nih.gov; do
    url_check "$u"
  done
  if [ "$PULL_TEST" = 1 ] && { has apptainer || has singularity; }; then
    local rt; has apptainer && rt=apptainer || rt=singularity
    local cache="${TMPDIR:-/tmp}/ff_inv_cache_$$"; mkdir -p "$cache"
    run 300 "apptainer pull test (compute node)" \
      "export APPTAINER_CACHEDIR=$cache SINGULARITY_CACHEDIR=$cache APPTAINER_TMPDIR=$cache SINGULARITY_TMPDIR=$cache; cd $cache && $rt pull -F alpine_test.sif docker://alpine:3.20 && $rt exec alpine_test.sif cat /etc/alpine-release && echo PULL+EXEC_OK"
    rm -rf "$cache"
  fi
  echo; echo "[compute probe finished $(date)]"
}

if [ "$MODE" = "compute" ]; then compute_probe; exit 0; fi

# ==============================================================================
# LOGIN-NODE SURVEY
# ==============================================================================
survey(){
  echo "fungiforge HPC inventory v$INVENTORY_VERSION"
  echo "generated: $(date)   by: $(id -un)   on: $(hostname -f 2>/dev/null || hostname)"
  echo "options: submit-test=$SUBMIT_TEST pull-test=$PULL_TEST partition='${PARTITION:-}' account='${ACCOUNT:-}'"
  echo "report: $OUT"

  # ---------------------------------------------------------------- 1. system
  hdr "1. SYSTEM"
  run 10 "os"                  "cat /etc/os-release 2>/dev/null | head -6; uname -srmo"
  run 10 "arch / glibc"        "uname -m; ldd --version 2>&1 | head -1; getconf LONG_BIT"
  run 10 "login node cpu/mem"  "nproc; free -g; uptime"
  run 10 "user / groups"       "id; echo umask=\$(umask); echo SHELL=\$SHELL; echo HOME=\$HOME"
  run 10 "sudo?"               "sudo -n true 2>&1 | head -2; echo rc=\$?"
  run 10 "motd (site rules)"   "head -60 /etc/motd 2>/dev/null"
  run 10 "locale / tmp"        "echo LANG=\$LANG LC_ALL=\$LC_ALL; echo TMPDIR=\${TMPDIR:-unset}; df -hP \${TMPDIR:-/tmp} 2>/dev/null"
  run 10 "ulimits (nextflow head process)" "ulimit -a"
  run 10 "cgroup limit on login session" "cat /sys/fs/cgroup/user.slice/user-\$(id -u).slice/memory.max 2>/dev/null; cat /sys/fs/cgroup/memory.max 2>/dev/null; systemctl show user-\$(id -u).slice -p MemoryMax -p CPUQuota 2>/dev/null"
  run 10 "long-running processes allowed on login? (screen/tmux/nohup present)" "command -v screen tmux nohup"

  # ---------------------------------------------------------------- 2. scheduler
  hdr "2. SCHEDULER"
  if has sinfo; then
    echo "scheduler: SLURM"
    run 20 "slurm version"         "sinfo --version; scontrol version 2>/dev/null"
    run 30 "cluster config"        "scontrol show config 2>/dev/null | grep -iE '^(ClusterName|SlurmctldHost|MaxJobCount|MaxArraySize|MaxSubmitJobs|DefMemPerCPU|DefMemPerNode|MaxMemPerNode|SchedulerType|SelectType|SelectTypeParameters|PriorityType|AccountingStorageType|JobSubmitPlugins|TmpFS|PrologFlags|ProctrackType|TaskPlugin|Licenses|EnforcePartLimits|MaxTasksPerNode)'"
    run 30 "partitions (name avail timelimit nodes cpus mem gres features)" "sinfo -o '%20P %6a %12l %6D %5c %8m %12G %40f' | sort -u"
    run 30 "partition details"     "scontrol show partition 2>/dev/null | grep -E 'PartitionName|AllowGroups|AllowAccounts|AllowQos|Default=|MaxTime|DefaultTime|MaxNodes|MaxCPUsPerNode|MaxMemPerNode|DefMemPerCPU|QoS=|State=|TotalCPUs|TotalNodes|OverSubscribe'"
    run 30 "does a 'short' / 'long' partition exist? (hpc_slurm.config default routing)" "sinfo -h -o '%P' | tr -d '*' | sort -u | grep -xE 'short|long' || echo 'NEITHER short nor long -> must pass --slurm_partition or edit conf/hpc_slurm.config'"
    run 30 "node inventory (per node type)" "sinfo -h -N -o '%c %m %G %f' | sort | uniq -c | sort -rn | head -40"
    run 30 "node states summary"   "sinfo -h -o '%P %t %D' | sort | uniq -c"
    run 30 "my associations (account/partition/qos/limits)" "sacctmgr -P show assoc user=\$(id -un) format=cluster,account,partition,qos,defaultqos,maxjobs,maxsubmit,maxtresperjob,grptres 2>&1 | head -30"
    run 30 "qos limits"            "sacctmgr -P show qos format=name,maxwall,maxtresperjob,maxtrespu,maxjobspu,maxsubmitpu,priority 2>&1 | head -40"
    run 20 "my queue right now"    "squeue -u \$(id -un) 2>&1 | head -10"
    run 20 "cluster load"          "squeue -h -o '%t' | sort | uniq -c"
    run 20 "default account/partition for me" "sacctmgr -Pn show user \$(id -un) format=defaultaccount 2>&1; sinfo -h -o '%P' | grep '\*' | tr -d '*' | sed 's/^/default partition: /'"
    run 10 "sbatch/srun/sacct present" "command -v sbatch srun sacct scancel"
  elif has qstat; then
    echo "scheduler: PBS/Torque/SGE (fungiforge ships hpc_slurm.config; a PBS profile will be needed)"
    run 30 "qstat -Q / queues"   "qstat -Q 2>&1 | head -40; qconf -sql 2>&1 | head -40"
    run 30 "pbsnodes summary"    "pbsnodes -a 2>/dev/null | grep -E 'resources_available.(ncpus|mem)|state' | sort | uniq -c | head -30"
  elif has bsub; then
    echo "scheduler: LSF (fungiforge ships hpc_slurm.config; an LSF profile will be needed)"
    run 30 "bqueues"             "bqueues 2>&1 | head -40"
    run 30 "lshosts"             "lshosts 2>&1 | head -40"
  else
    echo "scheduler: NONE DETECTED on PATH (sinfo/qstat/bsub missing) — check if a module must be loaded first"
  fi

  # ---------------------------------------------------------------- 3. modules
  hdr "3. MODULE SYSTEM"
  if init_modules; then
    run 10 "module system"       "echo LMOD_VERSION=\${LMOD_VERSION:-unset} MODULES_VERSION=\${MODULE_VERSION:-unset} LMOD_PKG=\${LMOD_PKG:-unset} MODULESHOME=\${MODULESHOME:-unset}"
    runm "module --version"      "module --version 2>&1 | head -3"
    run 10 "MODULEPATH"          "echo \"\$MODULEPATH\" | tr ':' '\n'"
    run 20 "modulepath dirs: exist / writable by me" "for d in \$(echo \"\$MODULEPATH\" | tr ':' ' '); do if [ -d \"\$d\" ]; then w=no; [ -w \"\$d\" ] && w=YES; echo \"\$d  exists  writable=\$w  owner=\$(stat -c '%U:%G' \"\$d\" 2>/dev/null)\"; else echo \"\$d  MISSING\"; fi; done"
    run 10 "site module roots (where apps are installed)" "echo \"\$MODULEPATH\" | tr ':' '\n' | xargs -I{} dirname {} 2>/dev/null | sort -u"
    runm "currently loaded"      "module list 2>&1"
    runm "modules: nextflow/java/containers/conda" "module --terse avail 2>&1 | grep -iE '^(nextflow|java|jdk|openjdk|temurin|corretto|apptainer|singularity|conda|anaconda|miniconda|miniforge|mamba|micromamba|python|aria2|git)' | head -60"
    runm "modules: bio tools fungiforge could reuse" "module --terse avail 2>&1 | grep -iE '^(busco|funannotate|antismash|interproscan|iprscan|genemark|repeatmasker|repeatmodeler|tetools|eggnog|kraken2|flye|medaka|dorado|blast|hmmer|diamond|augustus|braker|compleasm|sourmash|skani|fastani|bwa|minimap2|samtools|kraken)' | head -60"
    runm "module spider nextflow (Lmod only)" "module spider nextflow 2>&1 | head -40"
    runm "module spider apptainer/singularity" "module spider apptainer 2>&1 | head -25; module spider singularity 2>&1 | head -25"
    runm "module spider java" "module spider java 2>&1 | head -30; module spider openjdk 2>&1 | head -20"
    runm "module avail fungiforge (existing install?)" "module --terse avail fungiforge 2>&1 | head"
    runm "module count total" "module --terse avail 2>&1 | grep -vc ':$'"
  else
    echo "NO module system detected (module command unavailable and no init script found)."
    echo "  -> fungiforge's Lmod modulefile (share/modulefiles/fungiforge/0.1.0.lua) can't be used as-is; a wrapper/profile.d script would be needed."
  fi

  # ---------------------------------------------------------------- 4. tools on PATH
  hdr "4. TOOLS ON PATH (login node, current environment)"
  run 10 "PATH (deduplicated)" "echo \"\$PATH\" | tr ':' '\n' | awk '!seen[\$0]++'"
  sub "core requirements"
  tool nextflow -v; tool java -version; tool apptainer --version; tool singularity --version
  sub "alternatives / helpers"
  tool docker --version; tool podman --version; tool conda --version; tool mamba --version; tool micromamba --version
  tool python3 --version; tool pip3 --version; tool git --version; tool aria2c --version; tool curl --version; tool wget --version
  tool gcc --version; tool make --version; tool rsync --version; tool tar --version; tool unzip -v; tool lua -v; tool screen --version; tool tmux -V
  sub "java details"
  run 20 "java" "echo JAVA_HOME=\${JAVA_HOME:-unset}; java -version 2>&1; readlink -f \$(command -v java) 2>/dev/null"
  run 10 "other JDKs on disk" "ls -d /usr/lib/jvm/* /opt/java* /opt/jdk* /usr/java/* 2>/dev/null"
  sub "nextflow details"
  run 60 "nextflow" "echo NXF_HOME=\${NXF_HOME:-unset} NXF_VER=\${NXF_VER:-unset} NXF_OPTS=\${NXF_OPTS:-unset} NXF_WORK=\${NXF_WORK:-unset} NXF_TEMP=\${NXF_TEMP:-unset} NXF_OFFLINE=\${NXF_OFFLINE:-unset} NXF_SINGULARITY_CACHEDIR=\${NXF_SINGULARITY_CACHEDIR:-unset} NXF_APPTAINER_CACHEDIR=\${NXF_APPTAINER_CACHEDIR:-unset} NXF_CONDA_CACHEDIR=\${NXF_CONDA_CACHEDIR:-unset}; readlink -f \$(command -v nextflow) 2>/dev/null; nextflow info 2>&1 | head -20"
  run 10 "~/.nextflow" "ls -la \$HOME/.nextflow 2>/dev/null | head; du -sh \$HOME/.nextflow 2>/dev/null; cat \$HOME/.nextflow/scm 2>/dev/null | head -5"
  run 10 "nextflow plugins/assets cached" "ls \$HOME/.nextflow/plugins \$HOME/.nextflow/assets 2>/dev/null | head -20"
  run 10 "nextflow installed elsewhere?" "ls -l /usr/local/bin/nextflow /opt/nextflow* /apps/nextflow* /software/nextflow* 2>/dev/null; find /opt /apps /software /usr/local /cm/shared -maxdepth 4 -name 'nextflow' -type f 2>/dev/null | head"

  # ---------------------------------------------------------------- 5. containers
  hdr "5. CONTAINER RUNTIME"
  local RT=""; has apptainer && RT=apptainer; [ -z "$RT" ] && has singularity && RT=singularity
  if [ -n "$RT" ]; then
    echo "runtime: $RT  ($(command -v $RT))"
    run 20 "version / build config" "$RT --version; $RT buildcfg 2>/dev/null | grep -iE 'CONFIGDIR|SESSIONDIR|LIBEXECDIR|PACKAGE_VERSION' | head"
    run 20 "site config (readable parts)" "for f in /etc/apptainer/apptainer.conf /etc/singularity/singularity.conf /usr/local/etc/apptainer/apptainer.conf /usr/local/etc/singularity/singularity.conf; do [ -r \$f ] && { echo \"== \$f\"; grep -E '^(allow setuid|allow pid ns|user bind control|enable overlay|enable underlay|mount hostfs|mount tmp|mount home|bind path|enable fusemount|max loop devices|sessiondir max size|allow container|limit container|always use nv|config passwd|config group|shared loop devices|allow net users|allow net groups|allow kernel squashfs|allow setuid-mount)' \$f; }; done"
    run 10 "starter-suid (setuid mode) present?" "ls -l \$(dirname \$(readlink -f \$(command -v $RT)))/../libexec/$RT/bin/starter-suid 2>/dev/null || ls -l /usr/libexec/$RT/bin/starter-suid /usr/local/libexec/$RT/bin/starter-suid 2>/dev/null || echo 'no starter-suid found (rootless / user-namespace mode)'"
    run 10 "user namespaces (needed if no setuid)" "echo max_user_namespaces=\$(cat /proc/sys/user/max_user_namespaces 2>/dev/null); echo unprivileged_userns_clone=\$(cat /proc/sys/kernel/unprivileged_userns_clone 2>/dev/null); grep -E \"^\$(id -un):\" /etc/subuid /etc/subgid 2>/dev/null || echo 'no subuid/subgid entry for me (fakeroot builds unavailable)'"
    run 10 "container env vars" "env | grep -E '^(APPTAINER|SINGULARITY)_' | sort"
    run 20 "remote endpoints / registry logins" "$RT remote list 2>&1 | head -10"
    run 10 "cache dirs" "echo APPTAINER_CACHEDIR=\${APPTAINER_CACHEDIR:-unset} SINGULARITY_CACHEDIR=\${SINGULARITY_CACHEDIR:-unset}; du -sh \$HOME/.apptainer \$HOME/.singularity 2>/dev/null; ls \$HOME/.apptainer/cache \$HOME/.singularity/cache 2>/dev/null | head"
    run 10 "can I build .sif here? (fakeroot/--remote/proot)" "$RT build --help 2>&1 | grep -E 'fakeroot|remote|userns|proot' | head; command -v proot"
    run 30 "existing .sif images on shared paths" "find /opt /apps /software /shared /cvmfs/singularity.galaxyproject.org /cm/shared /usr/local -maxdepth 4 \( -name '*.sif' -o -name '*.simg' \) 2>/dev/null | head -30"
    if [ "$PULL_TEST" = 1 ]; then
      local cache="${TMPDIR:-/tmp}/ff_inv_cache_$$"; mkdir -p "$cache"
      run 300 "PULL TEST: docker://alpine:3.20 (login node)" \
        "export APPTAINER_CACHEDIR=$cache SINGULARITY_CACHEDIR=$cache APPTAINER_TMPDIR=$cache SINGULARITY_TMPDIR=$cache; cd $cache && $RT pull -F alpine_test.sif docker://alpine:3.20 && $RT exec alpine_test.sif cat /etc/alpine-release && $RT exec -B \$HOME alpine_test.sif ls -d \$HOME && echo PULL+EXEC+BIND_OK"
      rm -rf "$cache"
    else
      echo; echo "(pull test skipped — rerun with --pull-test to prove Docker Hub pulls work)"
    fi
  else
    echo "NO apptainer/singularity on PATH. Check section 3 for a module, or it may be absent entirely."
  fi
  sub "docker / podman (usually unavailable on HPC)"
  run 10 "docker" "docker info 2>&1 | head -5; groups | tr ' ' '\n' | grep -x docker"
  run 10 "podman" "podman info 2>&1 | head -5"

  # ---------------------------------------------------------------- 6. storage
  hdr "6. STORAGE"
  run 20 "all filesystems (type mount)" "mount | awk '{print \$5, \$3}' | grep -Ei 'lustre|gpfs|nfs|beegfs|ceph|xfs|ext4|zfs|panfs|weka|vast' | sort | uniq -c | sort -rn | head -40"
  run 20 "free space: home + candidate roots" "df -hP \$HOME /tmp /scratch /scratch/\$(id -un) /work /project /projects /data /shared /apps /software /opt /usr/local /cm/shared /sw /storage /nfs /cvmfs /lustre /gpfs /home 2>/dev/null | sort -u"
  run 20 "space in each root the DB (150-250 GB) + images (~30 GB) + work (100s GB) could live in" "df -hP 2>/dev/null | awk 'NR==1 || \$4 ~ /[0-9]T/ || (\$4 ~ /G/ && \$4+0 > 300)' | sort -k4 -h | tail -30"
  run 20 "quotas (user)" "quota -s -u \$(id -un) 2>&1 | head -20"
  run 40 "lustre quotas" "for m in \$(mount -t lustre 2>/dev/null | awk '{print \$3}'); do echo \"== \$m\"; lfs quota -h -u \$(id -un) \$m 2>&1 | head -5; for g in \$(id -Gn); do lfs quota -h -g \$g \$m 2>/dev/null | tail -2 | sed \"s/^/  group \$g: /\"; done; done"
  run 40 "gpfs quotas" "command -v mmlsquota >/dev/null && mmlsquota --block-size=auto 2>&1 | head -20"
  run 40 "beegfs quotas" "command -v beegfs-ctl >/dev/null && beegfs-ctl --getquota --uid \$(id -un) 2>&1 | head"
  run 10 "scratch policy hints" "ls -ld /scratch /scratch/\$(id -un) /work /project* /data 2>/dev/null; cat /scratch/README* /scratch/*POLICY* 2>/dev/null | head -20"
  run 30 "candidate SHARED INSTALL ROOTS — writable by me / by my groups" "for d in /opt /opt/apps /opt/software /apps /software /shared /shared/apps /usr/local /cm/shared /cm/shared/apps /sw /projects /project /work /data /storage /nfs /home/shared /group /groups /lab /labs \$(echo \"\$MODULEPATH\" | tr ':' '\n' | xargs -I{} dirname {} 2>/dev/null); do [ -d \"\$d\" ] || continue; w=no; [ -w \"\$d\" ] && w=YES; echo \"\$(stat -c '%A %U:%G' \"\$d\" 2>/dev/null)  writable=\$w  \$d\"; done | sort -u"
  run 30 "group-owned dirs I can write to under project/work/scratch roots" "for r in /project /projects /work /scratch /data /group /groups /lab /labs /storage /nfs /shared; do [ -d \$r ] && find \$r -maxdepth 2 -type d -writable 2>/dev/null | head -8; done"
  run 10 "paths containing SPACES (breaks funannotate/antiSMASH)" "for p in \$HOME \${TMPDIR:-/tmp} \$(df -P 2>/dev/null | awk 'NR>1{print \$6}'); do case \"\$p\" in *' '*) echo \"SPACE IN PATH: \$p\";; esac; done; echo '(none listed above = OK)'"
  run 10 "home inode usage" "df -iP \$HOME 2>/dev/null"

  # ---------------------------------------------------------------- 7. network
  hdr "7. OUTBOUND NETWORK (login node)"
  run 10 "proxy env" "env | grep -iE '^(http_proxy|https_proxy|ftp_proxy|no_proxy|HTTP_PROXY|HTTPS_PROXY|NO_PROXY)=' || echo '(no proxy vars)'"
  run 10 "dns" "cat /etc/resolv.conf 2>/dev/null | grep -vE '^#' | head -5; getent hosts github.com registry-1.docker.io 2>&1"
  sub "HTTP reachability (code time ip) — every host fungiforge fetches from"
  echo "  legend: any 3-digit code (200/301/401/403/404) = host REACHABLE; 000 = blocked/unreachable/timeout"
  for u in https://github.com/aduton1000/fungiforge \
           https://registry-1.docker.io/v2/ https://auth.docker.io/token https://production.cloudflare.docker.com \
           https://quay.io https://ghcr.io https://depot.galaxyproject.org/singularity/ \
           https://get.nextflow.io https://github.com/nextflow-io/nextflow/releases https://repo1.maven.org/maven2/ \
           https://conda.anaconda.org/bioconda https://conda.anaconda.org/conda-forge https://pypi.org/simple/ \
           https://busco-data.ezlab.org/v5/data/ https://dl.secondarymetabolites.org/releases/ https://unite.ut.ee \
           https://api.plutof.ut.ee https://files.plutof.ut.ee https://genome-idx.s3.amazonaws.com https://ftp.ncbi.nlm.nih.gov \
           https://osf.io http://eggnog5.embl.de/ https://ftp.ebi.ac.uk/pub/databases/interpro/ \
           https://www.dfam.org https://hgdownload.soe.ucsc.edu https://rvdb.dbi.udel.edu https://mardy.dide.ic.ac.uk; do
    url_check "$u"
  done
  run 30 "git clone possible?" "git ls-remote --heads https://github.com/aduton1000/fungiforge.git 2>&1 | head -3"
  run 20 "ssh to github (alt if https blocked)" "ssh -o BatchMode=yes -o ConnectTimeout=8 -T git@github.com 2>&1 | head -2"
  echo; echo "NOTE: compute nodes often have DIFFERENT (or no) internet access — see the compute-node probe (--submit-test)."

  # ---------------------------------------------------------------- 8. existing databases
  hdr "8. EXISTING REFERENCE DATABASES / LICENSES (reuse instead of re-downloading?)"
  run 10 "db-related env vars" "env | grep -iE '(BUSCO|FUNANNOTATE|ANTISMASH|EGGNOG|KRAKEN|INTERPRO|GENEMARK|AUGUSTUS|BLASTDB|REPEATMASKER|DFAM|FUNGIFORGE|NXF_)' | sort"
  run 120 "db dirs on shared roots (depth<=4)" "find /db /data /databases /shared /apps /software /opt /cvmfs /projects /project /work /scratch /ref /refs /references /bio /biodb /cm/shared \${EXTRA_DB_ROOT:-/nonexistent} -maxdepth 4 -type d \( -iname '*busco*' -o -iname '*funannotate*' -o -iname '*antismash*' -o -iname '*eggnog*' -o -iname '*kraken*' -o -iname '*unite*' -o -iname '*interpro*' -o -iname '*genemark*' -o -iname '*repeatmasker*' -o -iname '*dfam*' -o -iname '*pfam*' -o -iname '*refseq*' -o -iname '*ncbi*' -o -iname '*fungi*' -o -iname '*rvdb*' -o -iname '*sourmash*' -o -iname '*fungiforge*' \) 2>/dev/null | head -60"
  run 10 "cvmfs (galaxy refdata / singularity)" "ls /cvmfs 2>/dev/null; ls /cvmfs/data.galaxyproject.org/byhand 2>/dev/null | head -20"
  run 10 "GeneMark license" "ls -l \$HOME/.gm_key /opt/genemark*/gm_key* 2>/dev/null || echo 'no gm_key (funannotate falls back to non-GeneMark predictors; a free academic key at exon.gatech.edu improves annotation)'"
  run 10 "conda envs present" "conda env list 2>/dev/null | head -20; ls \$HOME/.conda/envs \$HOME/miniconda3/envs \$HOME/miniforge3/envs 2>/dev/null | head"
  run 10 ".condarc / channels" "cat \$HOME/.condarc /etc/conda/.condarc 2>/dev/null"

  # ---------------------------------------------------------------- 9. existing fungiforge
  hdr "9. EXISTING FUNGIFORGE INSTALL?"
  run 10 "env / files" "env | grep '^FUNGIFORGE' ; command -v fungiforge fungiforge-run; ls -l \$HOME/.fungiforge-env.sh 2>/dev/null"
  run 60 "checkouts" "find \$HOME /opt /apps /software /shared /projects /project /work /scratch -maxdepth 4 -type d -name 'fungiforge' 2>/dev/null | head"
  if [ -f "$(dirname "$SELF")/../main.nf" ] && has nextflow; then
    run 120 "nextflow config resolution from this checkout (hpc_slurm,apptainer)" "cd $(dirname "$SELF")/.. && nextflow config -profile hpc_slurm,apptainer 2>&1 | head -80"
  fi

  # ---------------------------------------------------------------- 10. test job
  hdr "10. COMPUTE-NODE TEST JOB"
  if [ "$SUBMIT_TEST" = 1 ] && has sbatch; then
    local probe_out="${OUT%.txt}_compute_probe.txt"
    local opts=(--job-name=ff_inventory --time=00:10:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=4G --output="$probe_out")
    [ -n "$PARTITION" ] && opts+=(--partition="$PARTITION")
    [ -n "$ACCOUNT" ]   && opts+=(--account="$ACCOUNT")
    local pflag=""; [ "$PULL_TEST" = 1 ] && pflag="--pull-test"
    echo "submitting: sbatch ${opts[*]} --wrap 'bash $SELF --compute-probe $pflag'"
    local sub; sub="$(sbatch "${opts[@]}" --wrap "bash '$SELF' --compute-probe $pflag" 2>&1)"; echo "$sub"
    local jid; jid="$(echo "$sub" | grep -oE '[0-9]+$')"
    if [ -n "$jid" ]; then
      echo "waiting up to ${WAIT_SECS}s for job $jid ..."
      local waited=0
      while [ $waited -lt "$WAIT_SECS" ]; do
        squeue -h -j "$jid" 2>/dev/null | grep -q . || break
        sleep 15; waited=$((waited+15))
        [ $((waited % 120)) -eq 0 ] && echo "  still waiting (${waited}s): $(squeue -h -j "$jid" -o '%T %r' 2>/dev/null)"
      done
      run 30 "job accounting" "sacct -j $jid -X -o JobID,State,ExitCode,Elapsed,Partition,NodeList,ReqMem,MaxRSS -P 2>&1"
      if [ -s "$probe_out" ]; then
        echo; echo "----- compute-node probe output ($probe_out) -----"; cat "$probe_out"; echo "----- end probe -----"
      else
        echo "probe output not present yet ($probe_out). Job still queued or failed to start."
        echo "When it finishes, append it:  cat '$probe_out' >> '$OUT'"
      fi
    fi
  else
    echo "skipped (rerun with --submit-test to survey a compute node: internet egress, apptainer, scratch, cgroups)"
  fi

  hdr "END OF REPORT"
  echo "Bring this file back: $OUT"
}

survey 2>&1 | tee "$OUT"
echo
echo "Report written: $OUT  ($(wc -l < "$OUT") lines)"
echo "Copy it back with:  scp $(id -un)@$(hostname -f 2>/dev/null || hostname):$OUT ."
