# fungiforge — HPC deployment

Native (no emulation) on a Linux SLURM cluster with Apptainer/Singularity.

## 1. Get the code + tools
```bash
git clone https://github.com/aduton1000/fungiforge && cd fungiforge
# needs: nextflow >=23.10, java 17+, apptainer (or singularity)
```

## 2. Build the base image (once)
```bash
apptainer build fungiforge.sif env/fungiforge.def
# the heavy public images pull automatically at run time (Funannotate, antiSMASH, tetools, …)
```

## 3. Fetch reference databases (once, hours)
```bash
export FUNGIFORGE_DB=/scratch/$USER/fungiforge_db
bash bin/fetch_references.sh          # aria2c: multi-connection, resumable
```

## 4. Site config + module
```bash
cp share/fungiforge-env.sh.example ~/.fungiforge-env.sh    # edit paths/partition/account
# expose the launcher via Lmod:
module use /path/to/fungiforge/share/modulefiles
module load fungiforge
```

## 5. Run
```bash
fungiforge-run --samplesheet samples.csv           # uses hpc_slurm,apptainer
# or directly:
nextflow run main.nf -profile hpc_slurm,singularity \
    --samplesheet samples.csv --data_dir "$FUNGIFORGE_DB" \
    --slurm_partition genomics --slurm_account pi_adu --run_interproscan true
```

## Notes
- Launch from a **writable** dir (Nextflow writes `.nextflow/`, `work/`); put `work/` on scratch.
- `--run_interproscan true` is the default-worth-enabling on HPC (too slow under Mac emulation).
- Partition routing is time-based by default (`short`/`long`); override with `--slurm_partition`.
- Bind mounts: `APPTAINER_BINDPATH` must include `$FUNGIFORGE_DB` and the work dir (see env example).
