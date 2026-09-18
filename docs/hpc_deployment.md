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
  images/antismash-ff-8.0.0-r3.sif   from env/antismash-ff.Dockerfile
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
sudo chown $USER:<group> /hpc/opt/fungiforge /hpc/data/fungiforge
sudo chmod 2775 /hpc/opt/fungiforge /hpc/data/fungiforge

git clone https://github.com/aduton1000/fungiforge /tmp/ff && cd /tmp/ff
bash bin/hpc_install.sh --prefix /hpc/opt/fungiforge --db /hpc/data/fungiforge \
     --group <group> --partition <partition> --profile-d --smoke
```

What it does: clones/updates `repo/`; `docker build` → `apptainer build … docker-daemon://`
for the two fungiforge images (native amd64, ~20 GB, 20–40 min); pre-pulls the public
images (flye, medaka, polypolish, busco, tetools, funannotate ≈ 15 GB) into
`images/cache` with Nextflow's file naming so runs never pull; creates `cli-env`;
writes `site.config`, `fungiforge-env.sh` and the `bin/` wrappers from the templates
in `share/`; sets group read/exec; writes `/etc/profile.d/fungiforge.sh` (with
`--profile-d`, via sudo) so every login gets the launchers on PATH; and with
`--smoke` runs the 15-stage DAG as `-stub-run` through SLURM.

The public images are pinned to versioned tags/digests in `conf/base.config`; the cache file
names follow Nextflow's convention (`staphb-flye-2.9.6.img`, `nextgenusfs-funannotate@sha256-….img`),
so a checkout that changes a pin needs the new file pulled (`fungiforge-fetch-refs images`, or
rerun the installer) — an older `…-latest.img` in the cache is simply not used. Set
`params.image_cache_dir` in `site.config` to the same directory as `apptainer.cacheDir` so
`provenance.json` can record the sha256 of every image a run used.

No Docker on the cluster? Build the two images on any Linux box with Docker
(`docker build --platform linux/amd64 …`, then `apptainer build x.sif docker-daemon://tag`),
copy the `.sif` files into `<prefix>/images/`, and rerun the installer with `--skip-images`.

## 2. Reference databases (one-time, hours, resumable)

```bash
/hpc/opt/fungiforge/bin/fungiforge-fetch-refs                 # every default step
/hpc/opt/fungiforge/bin/fungiforge-fetch-refs eggnog genomad  # selected steps
```

Each database is guarded by a `.done` marker, one failing step never aborts the others, and the
manifest (`$FUNGIFORGE_DB/MANIFEST.tsv`) records what succeeded. The two databases that are built
*through* their tool container (funannotate, antiSMASH) run under **Apptainer/Singularity** when
present and Docker otherwise, so no Docker daemon is needed on a cluster.

| Step | What it stages | Size |
|:--|:--|:--|
| `images` | every public image into the shared Apptainer cache | ~40 GB |
| `funannotate`, `antismash` | annotation and BGC databases (through their containers) | ~40 GB, ~9 GB |
| `eggnog` | eggNOG 5.0.2 (`eggnog.db`, `eggnog_proteins.dmnd`, taxa) — direct download | ~12 GB |
| `busco` | `fungi_odb10` plus the order/class lineages of `busco_lineages.tsv`, direct from the BUSCO data server | ~1 GB |
| `unite`, `kraken2`, `refseq_fungi` | ITS reference, contamination DB, sourmash signatures | ~15 GB |
| `fungamr` | resistance catalogue + reference proteins + derived panel | small |
| `markers`, `mlst` | type-material sets for CaM/BenA/TEF1/RPB2/LSU; PubMLST fungal schemes | ~100 MB |
| `dbcan`, `phibase`, `effectorp` | CAZyme HMMs, virulence proteins, EffectorP 3 | ~300 MB |
| `genomad` | geNomad database v1.9 | 0.8 GB |
| `benchmarks` | A1163 and NRRL 3357 reference genomes | ~80 MB |
| `genomes` *(on request)* | one reference genome per species of `novelty_genera.txt`, for genome-ANI novelty | tens of GB |
| `interproscan` *(on request)* | InterProScan data release matching the image tag (`--run_interproscan true` to use it) | 6.9 GB |

The multi-GB databases (eggNOG, Kraken 2, geNomad, UNITE) are fetched with `aria2c` when it is
installed, otherwise `wget -c`, otherwise `curl`. Install `aria2c` if you can: it is the only one of
the three that downloads in parallel chunks, and on a link that drops connections `curl` restarts a
transfer from byte 0 each retry, so the file grows and shrinks without ever finishing. Every archive
is integrity-checked before decompression and deleted if it fails, so a rerun starts clean. The step
only marks itself done when the decompressed files meet their expected sizes, and one fetch per data
directory runs at a time: a second one refuses to start rather than writing over the first, and
names the process holding the lock so you can stop it. Pass `--wait` to queue behind a running
fetch instead, which is how you line up a second database while a multi-hour download finishes.

Two resources are **licensed** and staged by hand (Appendix A of the upgrade plan):

- **GeneMark-ES** — unpack the academic tarball and pass `--genemark_dir <dir> --genemark_key <key>`
  (both the directory and the key's directory are bound into the container; the stage copies the
  key to a per-task `$HOME/.gm_key`, and prediction falls back to Augustus alone, recording why, if
  the key is unusable)
  to the run; the profiles bind the directory into the annotation container.
- **SignalP 6** — `bash bin/hpc_install.sh … --signalp <signalp-6.0*.fast.tar.gz>` builds a
  site-only image from it and points the extras stage at it in `site.config`.

The `fungiforge` CLI finds the pipeline through `FUNGIFORGE_HOME` (set by `fungiforge-env.sh`),
falling back to `$FUNGIFORGE_ROOT/repo`, then to its own source tree. It is pip-installed editable
from `repo/`, so `git -C <install_root>/repo pull` updates the CLI along with the pipeline; an
installation made before this was the case needs one `pip install -e` to catch up:

```bash
<install_root>/cli-env/bin/pip install -q --no-cache-dir -e <install_root>/repo
```
 Outside a login that sources
the site environment, set `FUNGIFORGE_HOME` to the checkout or the helper subcommands (`check`,
`validate`, `check-master`) cannot find `bin/`.

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
bash /hpc/opt/fungiforge/repo/bin/hpc_install.sh --prefix /hpc/opt/fungiforge --db /hpc/data/fungiforge --group <group> --ref v0.1.1
```
Existing `site.config` / `fungiforge-env.sh` are kept; images are rebuilt only with
`--rebuild-images` (or when the version in `nextflow.config` changes). Add `--signalp /path/to/signalp-6.0i.fast.tar.gz` to build the site-only SignalP 6 image for the extras stage (W2.6); it is written next to the other `.sif` files and wired into `site.config`.

## Notes
- Launch from a **writable** dir (Nextflow writes `.nextflow/`, `work/`, the log to CWD).
- DB and install paths must be space-free (funannotate/antiSMASH pass paths unquoted).
- `--run_interproscan true` and `--genome_id true` are on in the site config (native speed).
- Apptainer auto-mounts `$HOME`, `/tmp`, CWD only — anything else referenced in place
  must be in `FUNGIFORGE_BIND` (default `/hpc`).
