// Stage 07b — eggNOG-mapper (W2.5). Orthology-based functional annotation of the predicted
// proteins (eggnog-mapper 2.1 against the eggNOG 5 data staged by `fetch_references.sh eggnog`,
// fungal taxonomic scope). Its annotations table is handed to funannotate annotate --eggnog
// (the funannotate image ships no emapper). Skipped (empty table, status skipped) when the
// database is missing or still being staged; --skip_eggnog removes the stage. A failure here is
// never fatal: annotation continues without the eggNOG table and the stage records why.
process EGGNOG {
  tag { meta.id }
  label 'eggnog'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy', pattern: "*.{json,annotations}"
  input:  tuple val(meta), path(proteins)
  output: tuple val(meta), path("${meta.id}.emapper.annotations"), emit: annotations
          tuple val(meta), path("${meta.id}.eggnog.json"),         emit: json
  script:
  def data = params.eggnog_db ?: (params.data_dir ? params.data_dir + '/eggnog' : '')
  """
  source ff_status.sh; ff_init eggnog "${meta.id}" ${meta.id}.eggnog.json
  ff_version emapper -- bash -c "emapper.py --version 2>&1 | tail -1"
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  # The database must be COMPLETE, not merely present: a staging run that is still downloading or
  # decompressing leaves short files, and emapper would die on them mid-run. Sizes mirror the
  # thresholds bin/fetch_references.sh marks the step done at.
  DB_OK=0
  if [ -n "${data}" ]; then
    SZ_DB=\$(stat -c%s "${data}/eggnog.db" 2>/dev/null || echo 0)
    SZ_DMND=\$(stat -c%s "${data}/eggnog_proteins.dmnd" 2>/dev/null || echo 0)
    if [ "\$SZ_DB" -ge 8000000000 ] && [ "\$SZ_DMND" -ge 4000000000 ]; then DB_OK=1; fi
  fi
  if [ "\$DB_OK" = 1 ]; then
    sed 's/\\*\$//' ${proteins} > query.faa
    ff_run emapper --optional -- emapper.py -i query.faa --itype proteins -o ${meta.id} --output_dir . --data_dir "${data}" \\
        --cpu ${task.cpus} --tax_scope ${params.eggnog_tax_scope} --temp_dir "\$TMPDIR" --override
    # An emapper that died part-way leaves a truncated table; hand annotation nothing rather than half.
    [ "\$FF_RC" = 0 ] && [ -s ${meta.id}.emapper.annotations ] || : > ${meta.id}.emapper.annotations
    N=\$(grep -vc '^#' ${meta.id}.emapper.annotations || true)
  else
    ff_skip emapper "eggNOG data incomplete or absent under --data_dir/eggnog (eggnog.db=\${SZ_DB:-0}, eggnog_proteins.dmnd=\${SZ_DMND:-0}; bin/fetch_references.sh eggnog)"
    : > ${meta.id}.emapper.annotations; N=0
  fi
  printf '{"sample":"%s","stage":"eggnog","n_annotated":%s,"data_dir":"%s"}\\n' "${meta.id}" "\${N:-0}" "${data}" > ${meta.id}.eggnog.json
  ff_finalize
  """
  stub:
  ": > ${meta.id}.emapper.annotations; echo '{\"sample\":\"${meta.id}\",\"stage\":\"eggnog\",\"status\":\"skipped\"}' > ${meta.id}.eggnog.json"
}
