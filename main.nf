#!/usr/bin/env nextflow
// FungiForge — per-isolate fungal genomics from ONT (+ optional Illumina hybrid).
// Layer 1 of the fungiforge design: reads -> assembly -> polish -> decontam ->
// repeat-mask -> annotate -> identify -> resistance -> mobile -> BGC -> novelty
// -> eukaryote-extras -> per-isolate report + master table.
//
//   nextflow run main.nf -profile local,docker \
//     --samplesheet samples.csv --data_dir "/path/to/fungiforge_db"
nextflow.enable.dsl = 2

include { FUNGIFORGE } from './subworkflows/fungiforge.nf'

// ── help ──────────────────────────────────────────────────────────────────
def help() {
  log.info """
  ┌───────────────────────────────────────────────────────────────────────┐
  │  FungiForge ${workflow.manifest.version} — fungal ONT/hybrid genomics                     │
  └───────────────────────────────────────────────────────────────────────┘
  Usage:
    nextflow run main.nf -profile <exec>,<pkg> --samplesheet s.csv --data_dir DIR

  Profiles (compose one of each group):
    execution : local | hpc_slurm
    packaging : conda | mamba | docker | singularity | apptainer
    test      : synthetic/subsampled fixture

  Required:
    --samplesheet  CSV: sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season
    --data_dir     reference-database root (see bin/fetch_references.sh)

  Common options (see nextflow.config for all):
    --outdir       results dir            (default: results)
    --hybrid       auto|on|off            (Illumina polishing; default auto)
    --assembler    flye|canu|raven        (default flye)
    --busco_lineage auto|<odb10 lineage>  (default auto)
    --run_interproscan  true|false        (default false; heavy under emulation)
    --skip_bgc / --skip_mge / --skip_novelty / --skip_extras
  """.stripIndent()
}

// ── input parsing ───────────────────────────────────────────────────────────
def parse_row(row) {
  def meta = [ id         : row.sample,
               compartment: (row.compartment ?: 'NA'),
               facility   : (row.facility ?: 'NA'),
               season     : (row.season ?: 'NA') ]
  def files = [ ont: file(row.ont_fastq),
                r1 : (row.illumina_r1?.trim()) ? file(row.illumina_r1) : file("${projectDir}/assets/NO_R1"),
                r2 : (row.illumina_r2?.trim()) ? file(row.illumina_r2) : file("${projectDir}/assets/NO_R2") ]
  return tuple(meta, files)
}

// ── dispatch ─────────────────────────────────────────────────────────────────
workflow {
  if (params.help) { help(); return }

  // fail-loud guards (forge convention)
  if (!params.samplesheet)
    error "Missing --samplesheet (columns: sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season)"
  def pkg_containerised = ['docker','singularity','apptainer'].any { workflow.profile.contains(it) }
  if (!params.data_dir && !workflow.stubRun)
    log.warn "--data_dir not set: identification / resistance / BGC / novelty stages need reference DBs (see bin/fetch_references.sh)."
  if (params.basecall && !params.dorado_model)
    error "--basecall set but --dorado_model missing."

  samples = Channel.fromPath(params.samplesheet).splitCsv(header: true).map { parse_row(it) }
  FUNGIFORGE(samples)
}
