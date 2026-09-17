#!/usr/bin/env nextflow
// FungiForge — per-isolate fungal genomics from ONT, Illumina, or ONT+Illumina hybrid.
// Layer 1 of the fungiforge design: reads -> assembly -> polish -> decontam ->
// repeat-mask -> annotate -> identify -> resistance -> mobile -> BGC -> novelty
// -> eukaryote-extras -> per-isolate report + master table.
// Each isolate is routed by its inputs: long-read assembly (Flye, + Medaka, +
// optional Illumina hybrid polish) when ONT is present, else short-read assembly
// (SPAdes) from Illumina alone. See meta.assembly_mode below.
//
//   nextflow run main.nf -profile local,docker \
//     --samplesheet samples.csv --data_dir "/path/to/fungiforge_db"
nextflow.enable.dsl = 2

include { FUNGIFORGE } from './subworkflows/fungiforge.nf'

// ── help ──────────────────────────────────────────────────────────────────
def help() {
  log.info """
  ┌───────────────────────────────────────────────────────────────────────┐
  │  FungiForge ${workflow.manifest.version} — fungal ONT / Illumina / hybrid genomics        │
  └───────────────────────────────────────────────────────────────────────┘
  Usage:
    nextflow run main.nf -profile <exec>,<pkg> --samplesheet s.csv --data_dir DIR

  Profiles (compose one of each group):
    execution : local | hpc_slurm
    packaging : conda | mamba | docker | singularity | apptainer
    test      : synthetic/subsampled fixture

  Required:
    --samplesheet  CSV: sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season
                   Per isolate provide EITHER ont_fastq (long-read / hybrid) OR
                   illumina_r1+illumina_r2 (short-read only) — or all three (hybrid).
    --data_dir     reference-database root (see bin/fetch_references.sh)

  Input modes (auto-detected per row -> meta.assembly_mode):
    ONT only        ont_fastq set, Illumina blank      -> Flye + Medaka (resistance provisional)
    Hybrid          ont_fastq + illumina_r1/r2 set     -> Flye + Medaka + Polypolish (high conf)
    Illumina only   illumina_r1/r2 set, ont_fastq blank-> SPAdes short-read assembly (high conf)

  Common options (see nextflow.config for all):
    --outdir       results dir            (default: results)
    --hybrid       auto|on|off            (Illumina polishing of ONT; default auto)
    --assembler    flye|canu|raven        (long-read; default flye)
    --sr_assembler spades|megahit         (Illumina-only; default spades)
    --busco_lineage auto|<odb10 lineage>  (default auto)
    --run_interproscan  true|false        (default false; heavy — enable on a cluster)
    --skip_bgc / --skip_mge / --skip_novelty / --skip_extras
    --force_all    run the fungal stages even for isolates the gate would stop (non-fungal / failed QC)
  """.stripIndent()
}

// ── input parsing ───────────────────────────────────────────────────────────
// Each isolate needs EITHER ONT (long-read / hybrid) OR paired Illumina
// (short-read only). assembly_mode routes it through the right assembler in
// subworkflows/fungiforge.nf: 'longread' (ONT present, ± Illumina hybrid polish)
// or 'shortread' (Illumina only).
// Samplesheet paths: absolute, or relative to the launch dir (as before); a relative path
// that does not exist there is resolved against the samplesheet's own directory, so a
// sheet shipped next to its data (or the bundled test fixture) works from any launch dir.
// NB Nextflow's file() absolutises relative strings against the launch dir, so test the
// raw string, not the Path.
def resolve_path(p) {
  def s = p.toString()
  if (s.startsWith('/')) return file(s)
  def cwdf = file(s)
  if (cwdf.exists()) return cwdf
  def near = file(params.samplesheet).parent.resolve(s)
  return near.exists() ? near : cwdf
}

def parse_row(row) {
  def has_ont  = row.ont_fastq?.trim()
  def has_r1   = row.illumina_r1?.trim()
  def has_r2   = row.illumina_r2?.trim()
  if ((has_r1 as boolean) != (has_r2 as boolean))
    error "Sample '${row.sample}': Illumina must be PAIRED — provide both illumina_r1 and illumina_r2 (got r1='${row.illumina_r1}', r2='${row.illumina_r2}')."
  if (!has_ont && !(has_r1 && has_r2))
    error "Sample '${row.sample}': needs ont_fastq OR both illumina_r1+illumina_r2 (got ont='${row.ont_fastq}', r1='${row.illumina_r1}', r2='${row.illumina_r2}')."
  def meta = [ id           : row.sample,
               assembly_mode: (has_ont ? 'longread' : 'shortread'),
               compartment  : (row.compartment ?: 'NA'),
               facility     : (row.facility ?: 'NA'),
               season       : (row.season ?: 'NA') ]
  def files = [ ont: has_ont ? resolve_path(row.ont_fastq)   : file("${projectDir}/assets/NO_ONT"),
                r1 : has_r1  ? resolve_path(row.illumina_r1) : file("${projectDir}/assets/NO_R1"),
                r2 : has_r2  ? resolve_path(row.illumina_r2) : file("${projectDir}/assets/NO_R2") ]
  // Fail at launch, not inside the first tool: Nextflow stages a missing input as a dangling
  // link, and the read-QC tool then dies with an unhelpful "failed to open file".
  [ont_fastq: files.ont, illumina_r1: files.r1, illumina_r2: files.r2].each { col, p ->
    if (!p.exists())
      error "Sample '${row.sample}': ${col} file not found: '${row[col]}' (looked in the launch directory and next to the samplesheet)."
  }
  return tuple(meta, files)
}

// ── run metadata for provenance.json ─────────────────────────────────────────
// Written once per run into the work dir and handed to the PROVENANCE process (W0.2).
// The git commit is read here (head node) because the containers need not carry git.
def git_state() {
  def st = [commit: null, dirty: null]
  try {
    def p = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
    def out = p.text.trim(); p.waitFor()
    if (p.exitValue() == 0 && out) st.commit = out
    def q = ['git', '-C', projectDir.toString(), 'status', '--porcelain', '--untracked-files=no'].execute()
    def qo = q.text; q.waitFor()
    if (q.exitValue() == 0) st.dirty = qo.trim() ? true : false
  } catch (Exception _ignored) { }
  return st
}

def run_info_json() {
  def git = git_state()
  def host = null
  try { host = java.net.InetAddress.getLocalHost().getHostName() } catch (Exception _ignored) { }
  def containers = workflow.container instanceof Map ? workflow.container : [all: workflow.container?.toString()]
  def info = [
    pipeline : [name: workflow.manifest.name, version: workflow.manifest.version,
                repository: workflow.repository, revision: workflow.revision,
                commit_id: workflow.commitId ?: git.commit, git_dirty: git.dirty,
                script_id: workflow.scriptId, project_dir: projectDir.toString()],
    nextflow : [version: nextflow.version.toString(), build: nextflow.build],
    run      : [session_id: workflow.sessionId.toString(), run_name: workflow.runName,
                start: workflow.start.toString(), command_line: workflow.commandLine,
                launch_dir: workflow.launchDir.toString(), work_dir: workflow.workDir.toString(),
                profile: workflow.profile, container_engine: workflow.containerEngine,
                stub_run: workflow.stubRun, resume: workflow.resume, user: workflow.userName, host: host],
    containers: containers,
    params   : params,
  ]
  return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
}

// ── dispatch ─────────────────────────────────────────────────────────────────
workflow {
  if (params.help) { help(); return }

  // fail-loud guards (forge convention)
  if (!params.samplesheet)
    error "Missing --samplesheet (columns: sample,ont_fastq,illumina_r1,illumina_r2,compartment,facility,season; per row give ONT and/or paired Illumina)"
  if (!params.data_dir && !workflow.stubRun)
    log.warn "--data_dir not set: the DB_CHECK stage will stop the run (pass --allow_missing_db to run without reference databases; identification / decontamination / annotation / resistance / BGC then record skipped)."
  if (params.basecall && !params.dorado_model)
    error "--basecall set but --dorado_model missing."

  samples  = channel.fromPath(params.samplesheet).splitCsv(header: true).map { row -> parse_row(row) }
  run_info = channel.of(run_info_json()).collectFile(name: 'run_info.json', newLine: true)
  FUNGIFORGE(samples, run_info)

  // provenance: image identities (head node — needs the image cache / docker daemon; W0.2).
  // `params` is not visible inside the handler closure once the run has finished, so the
  // values it needs are captured here.
  def prov_out   = "${params.outdir}/pipeline_info/provenance.json"
  def prov_hash  = params.provenance_hash_images as boolean
  def prov_cache = params.image_cache_dir ?: System.getenv('NXF_APPTAINER_CACHEDIR') ?: System.getenv('NXF_SINGULARITY_CACHEDIR') ?: ''
  def prov_bin   = "${projectDir}/bin/make_provenance.py"
  workflow.onComplete = {
    def prov = file(prov_out)
    if (!workflow.success || !prov.exists()) return
    def cmd = ['python3', prov_bin, 'images', '--provenance', prov.toString(),
               '--engine', (workflow.containerEngine ?: 'none'), '--cache-dir', prov_cache, '--work-dir', workflow.workDir.toString()]
    if (!prov_hash) cmd << '--no-hash'
    try {
      def p = cmd.execute(); def out = p.text; p.waitFor()
      out.trim().eachLine { line -> log.info line }
      if (p.exitValue() != 0) log.warn "provenance image step exited ${p.exitValue()} (provenance.json kept without image identities)"
    } catch (Exception e) {
      log.warn "provenance image step could not run: ${e.message}"
    }
  }
  }
}
