// Stage 05b — the QC / verdict gate (W2.1, extended by W2.2). Runs ONLY for isolates the gate
// stops. Two gates use it: after read triage (reason `read_triage:<verdict>`: nothing from
// assembly onwards runs) and after assembly QC (reason `verdict:<verdict>` or `qc_pass:false`:
// repeat masking, annotation, identification, species-aware BUSCO, resistance, mobile, BGC,
// novelty and extras are not run). The isolate gets a `gate` stage JSON with status `skipped` and the reason, so its
// master row still appears (reason in the `gate` column, `gate:skipped` in stages_failed).
// --force_all disables both gates. The decision is made in the subworkflow from the stage
// JSONs; this process only records it.
process GATE {
  tag { meta.id }
  label 'report'
  publishDir { "${params.outdir}/${meta.id}/05_assembly_qc" }, mode: 'copy'
  input:  tuple val(meta), val(reason)
  output: tuple val(meta), path("${meta.id}.gate.json"), emit: json
  script:
  def skipped = reason.startsWith('read_triage')
    ? '["kmer", "assemble", "purge", "medaka", "polish", "decontam", "organelle", "assembly_qc", "repeatmask", "predict", "eggnog", "interproscan", "annotate", "identify", "busco_lineage", "resistance", "mobile", "bgc", "novelty", "extras"]'
    : '["repeatmask", "predict", "eggnog", "interproscan", "annotate", "identify", "busco_lineage", "resistance", "mobile", "bgc", "novelty", "extras"]'
  """
  python3 - <<'PY'
  import json
  json.dump({"sample": "${meta.id}", "stage": "gate", "status": "skipped", "reason": "${reason}",
             "skipped_stages": ${skipped},
             "tools": {}, "skipped_tools": {}, "versions": {},
             "note": "gate stopped this isolate (${reason}); the listed stages were not run (--force_all overrides)"},
            open("${meta.id}.gate.json", "w"), indent=2)
  PY
  """
  stub:
  def skipped = reason.startsWith('read_triage')
    ? '["kmer","assemble","purge","medaka","polish","decontam","organelle","assembly_qc","repeatmask","predict","eggnog","interproscan","annotate","identify","busco_lineage","resistance","mobile","bgc","novelty","extras"]'
    : '["repeatmask","predict","eggnog","interproscan","annotate","identify","busco_lineage","resistance","mobile","bgc","novelty","extras"]'
  """
  echo '{"sample":"${meta.id}","stage":"gate","status":"skipped","reason":"${reason}","skipped_stages":${skipped},"tools":{},"skipped_tools":{},"versions":{}}' > ${meta.id}.gate.json
  """
}
