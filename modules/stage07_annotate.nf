// Stage 07 — eukaryotic gene prediction + functional annotation (Funannotate).
// predict (GeneMark/Augustus/SNAP/GlimmerHMM via EVM) -> annotate (eggNOG, Pfam,
// dbCAN, MEROPS, antiSMASH handoff; InterProScan optional). Emits proteins + GBK.
process ANNOTATE {
  tag { meta.id }
  label 'annotate'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy'
  input:  tuple val(meta), path(masked)
  output: tuple val(meta), path("${meta.id}.proteins.faa"), emit: proteins
          tuple val(meta), path("${meta.id}.gbk"),          emit: gbk
          tuple val(meta), path("${meta.id}.annotate.json"),emit: json
  script:
  def ips = params.run_interproscan ? '' : '--no-interpro'
  """
  export FUNANNOTATE_DB="${params.data_dir}/funannotate"
  ${ params.genemark_key ? "cp ${params.genemark_key} ~/.gm_key || true" : "" }
  funannotate predict -i ${masked} -o fun -s "${meta.id}" \\
      --cpus ${task.cpus} --busco_db fungi ${ params.busco_lineage=='auto' ? '' : "--busco_seed_species ${params.busco_lineage}" } || true
  EGG=""; [ -n "\$(ls -A ${params.data_dir}/eggnog 2>/dev/null)" ] && EGG="--eggnog ${params.data_dir}/eggnog"
  funannotate annotate -i fun --cpus ${task.cpus} ${ips} \$EGG || true
  cp fun/predict_results/*.proteins.fa ${meta.id}.proteins.faa 2>/dev/null || \\
     cp fun/annotate_results/*.proteins.fa ${meta.id}.proteins.faa 2>/dev/null || touch ${meta.id}.proteins.faa
  cp fun/annotate_results/*.gbk ${meta.id}.gbk 2>/dev/null || cp fun/predict_results/*.gbk ${meta.id}.gbk 2>/dev/null || touch ${meta.id}.gbk
  printf '{"sample":"%s","stage":"annotate","proteins":"%s.proteins.faa","interproscan":%s}\\n' \\
    "${meta.id}" "${meta.id}" "${params.run_interproscan.toString()}" > ${meta.id}.annotate.json
  """
  stub:
  "touch ${meta.id}.proteins.faa ${meta.id}.gbk ${meta.id}.annotate.json"
}
