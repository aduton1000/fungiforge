# fungiforge — HPC deployment (SLURM + Apptainer, site-wide)

Native linux/amd64 on a SLURM cluster with Apptainer/Singularity. One central,
read-only-for-users install that every user invokes as `fungiforge-run`. Works
**without a module system** (site env sourced from `/etc/profile.d`); an Lmod
modulefile is also shipped for sites that have one.

Layout (mirrors callforge / captureforge on the same cluster):

```
<prefix>/                      e.g. /hpc/opt/fungiforge
  repo/                        git checkout (pinned --ref)
  images/fungiforge-0.1.0.sif  built here from env/Dockerfile (NOT on Docker Hub)
  images/antismash-ff-8.0.0-r2.sif   from env/antismash-ff.Dockerfile
  images/cache/                shared Nextflow image cache (public images pre-pulled)
  cli-env/                     conda env with the `fungiforge` Python CLI
  site.config                  nextflow -c overrides: caps, .sif paths, partition
  fungiforge-env.sh            site env (sourced by wrappers + /etc/profile.d)
  bin/fungiforge, fungiforge-run, fungiforge-preflight, fungiforge-fetch-refs
<db>/                          e.g. /hpc/data/fungiforge — reference DBs (150–250 GB)
```

## 0. Survey the cluster first (read-only)

```bash
bash bin/hpc_inventory.sh --submit-test --pull-test     # one report file to review
bash bin/hpc_forge_audit.sh                             # state of existing *forge installs
```

Requirements the survey checks: Nextflow ≥ 23.10 + Java 17+, Apptainer/Singularity
able to pull from Docker Hub (login **and** compute nodes), Docker (to build the two
fungiforge images; or build elsewhere and copy the `.sif`), SLURM, ≥ 300 GB on a
shared filesystem whose path contains **no spaces**, outbound HTTPS to the hosts in
`bin/fetch_references.sh`.

## 1. Install (idempotent — rerun to update)

```bash
# one-time, if the roots are root-owned:
sudo mkdir -p /hpc/opt/fungiforge /hpc/data/fungiforge
sudo chown $USER:hpcusers /hpc/opt/fungiforge /hpc/data/fungiforge
sudo chmod 2775 /hpc/opt/fungiforge /hpc/data/fungiforge

git clone https://github.com/aduton1000/fungiforge /tmp/ff && cd /tmp/ff
bash bin/hpc_install.sh --prefix /hpc/opt/fungiforge --db /hpc/data/fungiforge \
     --group hpcusers --partition global --profile-d --smoke
```

What it does: clones/updates `repo/`; `docker build` → `apptainer build … docker-daemon://`
for the two fungiforge images (native amd64, ~20 GB, 20–40 min); pre-pulls the public
images (flye, medaka, polypolish, busco, tetools, funannotate ≈ 15 GB) into
`images/cache` with Nextflow's file naming so runs never pull; creates `cli-env`;
writes `site.config`, `fungiforge-env.sh` and the `bin/` wrappers from the templates
in `share/`; sets group read/exec; writes `/etc/profile.d/fungiforge.sh` (with
`--profile-d`, via sudo) so every login gets the launchers on PATH; and with
`--smoke` runs the 15-stage DAG as `-stub-run` through SLURM.

No Docker on the cluster? Build the two images on any Linux box with Docker
(`docker build --platform linux/amd64 …`, then `apptainer build x.sif docker-daemon://tag`),
copy the `.sif` files into `<prefix>/images/`, and rerun the installer with `--skip-images`.

## 2. Reference databases (one-time, hours, resumable)

```bash
fungiforge-fetch-refs                       # all steps into $FUNGIFORGE_DB
fungiforge-fetch-refs unite fungamr busco   # or selected steps
```

Some steps run downloader tools inside containers via `docker run` (antiSMASH,
funannotate, eggNOG) — on a cluster where Docker is available to the installer this
works as-is; files land root-owned but world-readable. Optional: a free academic
GeneMark key in `~/.gm_key` (funannotate uses it if present).

## 3. Run (any user)

```bash
mkdir -p ~/runs/batch1 && cd ~/runs/batch1               # run state lives in CWD
fungiforge samplesheet --ont 'reads/*.fastq.gz' --compartment AIR --facility SITE_A --season WET -o samples.csv
fungiforge-preflight --input reads/ --outdir preflight    # ONT triage before a real run
fungiforge-run --samplesheet samples.csv --outdir results  # -resume to continue
```

`fungiforge-run` composes `-profile hpc_slurm,apptainer`, adds `-c site.config`,
`--data_dir $FUNGIFORGE_DB`, the partition/account from the env, binds the DB +
work dir (+ `FUNGIFORGE_BIND`) into task containers, points Nextflow at the shared
image cache, and refuses to start from a non-writable directory. Pass `-profile …`
yourself to override (e.g. `-profile test,hpc_slurm,apptainer -stub-run`).

Keep the head process alive across logouts: run inside `tmux`, or submit it as a
job: `sbatch --wrap 'fungiforge-run --samplesheet samples.csv' --time=7-0 --cpus-per-task=2 --mem=8G`.

## 4. Site tuning

- `site.config` (applied last, wins): `params.max_cpus/max_memory/max_time`, per-label
  `cpus/memory/time` for the big stages, `executor.queueSize`, `apptainer.pullTimeout`,
  `params.fungiforge_image/antismash_image`, `params.slurm_partition`.
- `fungiforge-env.sh`: partition/account, `FUNGIFORGE_WORK` (scratch), `FUNGIFORGE_BIND`.
- `conf/hpc_slurm.config`: partition is **unset by default** → cluster default partition;
  set `--slurm_partition` / `FUNGIFORGE_PARTITION` where the default is wrong.
- Lmod sites: `module use <prefix>/repo/share/modulefiles && module load fungiforge`.

## 5. Update

```bash
bash /hpc/opt/fungiforge/repo/bin/hpc_install.sh --prefix /hpc/opt/fungiforge --db /hpc/data/fungiforge --group hpcusers --ref v0.1.1
```
Existing `site.config` / `fungiforge-env.sh` are kept; images are rebuilt only with
`--rebuild-images` (or when the version in `nextflow.config` changes).

## Notes
- Launch from a **writable** dir (Nextflow writes `.nextflow/`, `work/`, the log to CWD).
- DB and install paths must be space-free (funannotate/antiSMASH pass paths unquoted).
- `--run_interproscan true` and `--genome_id true` are on in the site config (native speed).
- Apptainer auto-mounts `$HOME`, `/tmp`, CWD only — anything else referenced in place
  must be in `FUNGIFORGE_BIND` (default `/hpc`).
