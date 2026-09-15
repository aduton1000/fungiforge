// Stage 15 (run-level) — provenance.json: pipeline commit, Nextflow version, effective
// parameters, process->container map, database manifest, and every sample's per-stage
// status / exit codes / tool versions, consolidated into one tool-version table.
// Runs once after every isolate's report. main.nf's workflow.onComplete then adds the
// image identities (sha256 of the .sif/.img files that ran, or docker repo digests) —
// that step needs the head node, where the image cache and docker daemon are visible.
process PROVENANCE {
  tag 'run'
  label 'report'
  publishDir { "${params.outdir}/pipeline_info" }, mode: 'copy'
  input:  path(run_info)
          path(db_manifest)
          path(jsons,   stageAs: 'stage_json/*')
          path(masters, stageAs: 'master/*')
  output: path('provenance.json'), emit: provenance
  script:
  """
  python3 ${projectDir}/bin/make_provenance.py aggregate --run-info ${run_info} --db-manifest ${db_manifest} \\
      --jsons stage_json/*.json --out provenance.json
  python3 ${projectDir}/bin/make_provenance.py validate provenance.json
  """
  stub:
  """
  python3 ${projectDir}/bin/make_provenance.py aggregate --run-info ${run_info} --db-manifest ${db_manifest} \\
      --jsons stage_json/*.json --out provenance.json
  """
}
