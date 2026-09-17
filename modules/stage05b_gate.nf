// Stage 05b — the QC / verdict gate (W2.1). Runs ONLY for isolates the gate stops: an
// isolate whose contigs are not fungal (Kraken2 verdict non_fungal or human) or whose
// assembly failed QC (qc_pass false) gets a `gate` stage JSON with status `skipped` and
// the reason, so its master row still appears (with the reason in the `gate` column and
// `gate:skipped` in stages_failed) while repeat masking, annotation, identification,
// resistance, mobile, BGC, novelty and extras are not run for it. --force_all disables
// the gate. The decision itself is made in the subworkflow from the decontam and
// assembly-QC JSONs; this process only records it.
process GATE {
  tag { meta.id }
  label 'report'
  publishDir { "${params.outdir}/${meta.id}/05_assembly_qc" }, mode: 'copy'
  input:  tuple val(meta), val(reason)
  output: tuple val(meta), path("${meta.id}.gate.json"), emit: json
  script:
  """
  python3 - <<'PY'
  import json
  json.dump({"sample": "${meta.id}", "stage": "gate", "status": "skipped", "reason": "${reason}",
             "skipped_stages": ["repeatmask", "annotate", "identify", "resistance", "mobile", "bgc", "novelty", "extras"],
             "tools": {}, "skipped_tools": {}, "versions": {},
             "note": "gate stopped this isolate after assembly QC (${reason}); downstream fungal stages not run (--force_all overrides)"},
            open("${meta.id}.gate.json", "w"), indent=2)
  PY
  """
  stub:
  """
  echo '{"sample":"${meta.id}","stage":"gate","status":"skipped","reason":"${reason}","skipped_stages":["repeatmask","annotate","identify","resistance","mobile","bgc","novelty","extras"],"tools":{},"skipped_tools":{},"versions":{}}' > ${meta.id}.gate.json
  """
}
