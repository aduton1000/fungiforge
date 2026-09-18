// Stage 07c — InterProScan (W2.5). Domain / family / GO annotation of the predicted proteins in
// the official InterProScan image with the data release staged by `fetch_references.sh
// interproscan` (bound at /opt/interproscan/data by the profile). The XML feeds funannotate
// annotate --iprscan. Heavy (hours per genome): off by default (--run_interproscan true);
// skipped with an empty XML when the data release is not staged.
process INTERPROSCAN {
  tag { meta.id }
  label 'interproscan'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy', pattern: "*.{json,xml}"
  input:  tuple val(meta), path(proteins)
  output: tuple val(meta), path("${meta.id}.iprscan.xml"),       emit: xml
          tuple val(meta), path("${meta.id}.interproscan.json"), emit: json
  script:
  """
  source ff_status.sh; ff_init interproscan "${meta.id}" ${meta.id}.interproscan.json
  # the official image keeps interproscan.sh off PATH (its ENTRYPOINT, cleared by the docker profile)
  IPS=\$(command -v interproscan.sh || echo /opt/interproscan/interproscan.sh)
  ff_version interproscan -- bash -c "\$IPS --version 2>&1 | grep -m1 -i version"
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  if [ -d /opt/interproscan/data/pfam ] || [ -d /opt/interproscan/data/panther ]; then
    sed 's/\\*\$//' ${proteins} > query.faa
    ff_run interproscan -- "\$IPS" -i query.faa -f XML -o ${meta.id}.iprscan.xml -dp -goterms -iprlookup -pa \\
        --cpu ${task.cpus} -T "\$TMPDIR" ${params.interproscan_apps ? "-appl ${params.interproscan_apps}" : ''}
    N=\$(grep -c '<protein>' ${meta.id}.iprscan.xml || true)
  else
    ff_skip interproscan "InterProScan data not mounted at /opt/interproscan/data — stage it with 'fetch_references.sh interproscan' (expected at <data_dir>/interproscan/data) or pass --interproscan_data"
    : > ${meta.id}.iprscan.xml; N=0
  fi
  printf '{"sample":"%s","stage":"interproscan","n_proteins_scanned":%s}\\n' "${meta.id}" "\${N:-0}" > ${meta.id}.interproscan.json
  ff_finalize
  """
  stub:
  ": > ${meta.id}.iprscan.xml; echo '{\"sample\":\"${meta.id}\",\"stage\":\"interproscan\",\"status\":\"skipped\"}' > ${meta.id}.interproscan.json"
}
