// Stage 17 (run-level) — cohort summary (W3.3). Merges every isolate's master row into one
// results/04_summary/master_fungi.tsv in schema order, validates it against
// fungiforge/resources/master_schema.json (the Layer-2 contract: append-only columns, declared
// types, NA for anything not measured), renders cohort_report.html (composition, gates,
// resistance, mycotoxins, metric ranges and the stage-16 clusters/tree/clonal groups/BGC
// families), writes cohort_summary.json, and exports MultiQC custom content — running MultiQC
// itself when it is in the image. Runs after every REPORT; a schema error fails the stage
// (loudly, at the end of the run) without touching the per-isolate results.
process SUMMARY {
  tag 'run'
  label 'report'
  publishDir { "${params.outdir}/04_summary" }, mode: 'copy'
  input:  path(masters, stageAs: 'master/*')
          path(cohort_json)
          path(provenance)
  output: path('master_fungi.tsv'),     emit: master
          path('cohort_report.html'),   emit: report
          path('cohort_summary.json'),  emit: summary
          path('multiqc*'),             emit: multiqc, optional: true
  script:
  def coh = cohort_json.name.startsWith('NO_') ? '' : "--cohort ${cohort_json}"
  def prov = provenance.name.startsWith('NO_') ? '' : "--provenance ${provenance}"
  """
  source ff_status.sh; ff_init summary run summary_status.json --best-effort
  ff_version python -- python3 --version
  ff_run cohort_report -- cohort_report.py --masters master --outdir . ${coh} ${prov} \\
      --schema "${params.master_schema}" --version ${workflow.manifest.version}
  if command -v multiqc >/dev/null 2>&1 && [ -d multiqc ]; then
    ff_run multiqc --optional -- bash -c "multiqc --force --quiet --config multiqc/multiqc_config.yaml --filename multiqc_report.html --outdir multiqc_report multiqc > multiqc.log 2>&1"
  else
    ff_skip multiqc "multiqc not in this image; the custom-content JSON under multiqc/ is written for an external run"
  fi
  ff_finalize
  """
  stub:
  """
  printf 'sample\\n' > master_fungi.tsv
  echo '<html>cohort</html>' > cohort_report.html
  echo '{"n_isolates":0}' > cohort_summary.json
  mkdir -p multiqc
  """
}
