// Stage 14 — per-isolate report + master table row + provenance. Merges every
// per-stage result.json into a self-contained HTML report and appends the
// isolate's row to results/04_summary/master_fungi.tsv (the Layer-2 handoff).
process REPORT {
  tag { meta.id }
  label 'report'
  publishDir { "${params.outdir}/${meta.id}/14_report" }, mode: 'copy', pattern: "*.html"
  publishDir { "${params.outdir}/04_summary" },            mode: 'copy', pattern: "*.tsv"
  input:  tuple val(meta), path(jsons)
  output: tuple val(meta), path("${meta.id}.report.html"), emit: report
          path("${meta.id}.master.tsv"),                   emit: master
  script:
  """
  python3 ${projectDir}/bin/make_report.py --sample "${meta.id}" \\
      --compartment "${meta.compartment}" --facility "${meta.facility}" --season "${meta.season}" \\
      --jsons ${jsons} \\
      --html ${meta.id}.report.html --master ${meta.id}.master.tsv \\
      || { echo -e "sample\\tcompartment\\tspecies" > ${meta.id}.master.tsv; echo -e "${meta.id}\\t${meta.compartment}\\tNA" >> ${meta.id}.master.tsv; echo "<html><body>${meta.id}</body></html>" > ${meta.id}.report.html; }
  """
  stub:
  "echo '<html>${meta.id}</html>' > ${meta.id}.report.html; printf 'sample\\n${meta.id}\\n' > ${meta.id}.master.tsv"
}
