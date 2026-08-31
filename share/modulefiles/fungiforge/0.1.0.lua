-- Lmod modulefile for fungiforge 0.1.0 — `module load fungiforge`
whatis("fungiforge 0.1.0 — fungal ONT/hybrid genomics pipeline")
help([[
Loads fungiforge onto PATH and sources ~/.fungiforge-env.sh.
Run:  fungiforge-run --samplesheet samples.csv
Needs: nextflow (>=23.10), apptainer/singularity, java 17+.
]])

local home = os.getenv("FUNGIFORGE_HOME") or "/opt/fungiforge"
prepend_path("PATH", pathJoin(home, "share/bin"))
setenv("FUNGIFORGE_HOME", home)

-- site config (paths, partition) is sourced by fungiforge-run from ~/.fungiforge-env.sh
depends_on("nextflow")
