// Stage 03a — ONT polishing with Medaka. Always runs. The optional short-read
// (Illumina) Polypolish hybrid pass is Stage 03b (stage03b_srpolish.nf), split
// out because Medaka and Polypolish live in different containers.
process MEDAKA {
  tag { meta.id }
  label 'polish'
  publishDir { "${params.outdir}/${meta.id}/03_polish" }, mode: 'copy'
  input:  tuple val(meta), path(assembly), path(ont)
  output: tuple val(meta), path("${meta.id}.medaka.fasta"), emit: assembly
          tuple val(meta), path("${meta.id}.medaka.json"),  emit: json
  script:
  """
  source ff_status.sh; ff_init medaka "${meta.id}" ${meta.id}.medaka.json
  ff_version medaka -- medaka --version
  mkdir -p medaka
  # A failed polish is a failed stage: silently carrying the unpolished assembly forward
  # would make every downstream resistance call untrustworthy without saying so.
  ff_run medaka -- medaka_consensus -i ${ont} -d ${assembly} -o medaka -t ${task.cpus} -m ${params.dorado_model}
  ff_run medaka_output -- test -s medaka/consensus.fasta
  cp medaka/consensus.fasta ${meta.id}.medaka.fasta
  printf '{"sample":"%s","stage":"medaka","model":"%s"}\\n' "${meta.id}" "${params.dorado_model}" > ${meta.id}.medaka.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.medaka.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"medaka\"}' > ${meta.id}.medaka.json"
}
