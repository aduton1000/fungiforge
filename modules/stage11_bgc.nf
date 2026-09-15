// Stage 11 — biosynthetic gene clusters (fungiSMASH / antiSMASH 8 fungal mode).
// NRPS/PKS/terpene/RiPP/hybrid clusters; downstream BiG-SCAPE networks + MIBiG
// (known mycotoxin/antibiotic vs novel) run at the comparative layer across isolates.
process BGC {
  tag { meta.id }
  label 'bgc'
  publishDir { "${params.outdir}/${meta.id}/11_bgc" }, mode: 'copy'
  input:  tuple val(meta), path(gbk)
  output: tuple val(meta), path("${meta.id}.bgc.json"),      emit: json
          tuple val(meta), path("${meta.id}.regions.gbk"),   emit: regions, optional: true
  script:
  """
  source "${projectDir}/bin/ff_status.sh"; ff_init bgc "${meta.id}" ${meta.id}.bgc.json --best-effort
  # Best-effort stage: an antiSMASH failure is RECORDED (status failed, n_clusters null -> NA
  # in master_fungi.tsv) but does not stop a multi-hour run.
  mkdir -p as; STATUS=skipped
  if [ "${params.skip_bgc}" = "true" ]; then ff_skip antismash "disabled (--skip_bgc)"
  elif ! grep -q "^LOCUS" ${gbk} 2>/dev/null || ! grep -qE "^     CDS |/translation=" ${gbk} 2>/dev/null; then
    ff_skip antismash "GBK is empty or has no annotated CDS"
  else
    ff_run antismash --optional -- bash -c "antismash --taxon fungi --output-dir as --genefinding-tool none --cpus ${task.cpus} --databases '${params.antismash_db ?: params.data_dir + '/antismash'}' ${gbk} > as.log 2>&1"
    if [ "\$FF_RC" -eq 0 ]; then STATUS=ok; else STATUS=failed; tail -20 as.log >&2; fi
  fi
  if ls as/*.region*.gbk >/dev/null 2>&1; then cat as/*.region*.gbk > ${meta.id}.regions.gbk; fi
  ff_run bgc_summary -- python3 ${projectDir}/bin/bgc_summary.py --sample "${meta.id}" --as-dir as --status "\$STATUS" \\
      --log as.log --out ${meta.id}.bgc.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"bgc\"}' > ${meta.id}.bgc.json"
}
