// Stage 08b — species-aware BUSCO (W2.3). Stage 05 scores completeness against fungi_odb10
// before the species is known (the QC gate); this stage re-scores the nuclear assembly with
// the lineage stage 08 chose from the species call (fungiforge/resources/busco_lineages.tsv,
// e.g. eurotiales_odb10 for Aspergillus). Skipped (status skipped, nulls) when the chosen
// lineage is the one stage 05 already used, when it is not staged under <data_dir>/busco,
// or when --busco_lineage names a fixed lineage. Emits <id>.busco_lineage.json.
process BUSCO_LINEAGE {
  tag { meta.id }
  label 'busco'
  publishDir { "${params.outdir}/${meta.id}/08_identify" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear), path(lineage_file)
  output: tuple val(meta), path("${meta.id}.busco_lineage.json"), emit: json
  script:
  def qc_lineage = params.busco_lineage == 'auto' ? 'fungi_odb10' : params.busco_lineage
  def data = params.data_dir ?: ''
  """
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  source ff_status.sh; ff_init busco_lineage "${meta.id}" ${meta.id}.busco_lineage.json
  ff_version busco -- busco --version
  ff_version compleasm -- compleasm --version
  LINEAGE=\$(head -1 ${lineage_file} | tr -d '[:space:]')
  STAGED=""
  for d in "${data}/busco/\$LINEAGE" "${data}/busco/lineages/\$LINEAGE" "${data}/busco/busco_downloads/lineages/\$LINEAGE"; do
    [ -d "\$d" ] && STAGED="\$d" && break
  done
  if [ "${params.busco_lineage}" != "auto" ]; then
    ff_skip busco "--busco_lineage ${params.busco_lineage} is fixed; stage 05 already scored it"
    busco_lineage.py --sample "${meta.id}" --lineage "\$LINEAGE" --qc-lineage "${qc_lineage}" --skipped fixed_lineage --out ${meta.id}.busco_lineage.json
  elif [ -z "\$LINEAGE" ] || [ "\$LINEAGE" = "${qc_lineage}" ]; then
    ff_skip busco "species-aware lineage is ${qc_lineage}, the lineage of stage 05"
    busco_lineage.py --sample "${meta.id}" --lineage "\$LINEAGE" --qc-lineage "${qc_lineage}" --skipped same_as_qc --out ${meta.id}.busco_lineage.json
  elif [ -z "\$STAGED" ]; then
    ff_skip busco "lineage \$LINEAGE not staged under --data_dir/busco (bin/fetch_references.sh busco)"
    busco_lineage.py --sample "${meta.id}" --lineage "\$LINEAGE" --qc-lineage "${qc_lineage}" --skipped not_staged --out ${meta.id}.busco_lineage.json
  else
    OK=0
    if command -v compleasm >/dev/null 2>&1; then
      ff_run compleasm --optional -- compleasm run -a ${nuclear} -o compleasm -l "\$LINEAGE" -L "${data}/busco" -t ${task.cpus}
      [ "\$FF_RC" -eq 0 ] && OK=1
    fi
    if [ "\$OK" -eq 0 ]; then
      ff_run busco -- busco -i ${nuclear} -o busco -l "\$LINEAGE" -m genome -c ${task.cpus} --download_path "${data}/busco" --offline
    fi
    busco_lineage.py --sample "${meta.id}" --lineage "\$LINEAGE" --qc-lineage "${qc_lineage}" --compleasm-dir compleasm --busco-dir busco --out ${meta.id}.busco_lineage.json
  fi
  ff_finalize
  """
  stub:
  """
  echo '{"sample":"${meta.id}","stage":"busco_lineage","lineage":"fungi_odb10","busco_complete":null,"status":"skipped"}' > ${meta.id}.busco_lineage.json
  """
}
