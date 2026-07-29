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
  if [ "${params.skip_bgc}" != "true" ]; then
    antismash --taxon fungi --output-dir as --genefinding-tool none --cpus ${task.cpus} \\
        --databases "${params.data_dir}/antismash" ${gbk} 2>as.log || true
    cat as/*.region*.gbk > ${meta.id}.regions.gbk 2>/dev/null || true
  fi
  python3 ${projectDir}/bin/bgc_summary.py --sample "${meta.id}" --as-dir as \\
      --out ${meta.id}.bgc.json || echo '{"sample":"${meta.id}","stage":"bgc","clusters":[]}' > ${meta.id}.bgc.json
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"bgc\"}' > ${meta.id}.bgc.json"
}
