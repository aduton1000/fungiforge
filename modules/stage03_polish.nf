// Stage 03a — ONT polishing with Medaka. Always runs. The optional short-read
// (Illumina) Polypolish hybrid pass is Stage 03b (stage03b_srpolish.nf), split
// out because Medaka and Polypolish live in different containers.
process MEDAKA {
  tag { meta.id }
  label 'polish'
  publishDir { "${params.outdir}/${meta.id}/03_polish" }, mode: 'copy'
  input:  tuple val(meta), path(assembly), path(ont)
  output: tuple val(meta), path("${meta.id}.medaka.fasta"), emit: assembly
  script:
  """
  mkdir -p medaka
  medaka_consensus -i ${ont} -d ${assembly} -o medaka -t ${task.cpus} -m ${params.dorado_model} || true
  if [ -s medaka/consensus.fasta ]; then cp medaka/consensus.fasta ${meta.id}.medaka.fasta
  else cp ${assembly} ${meta.id}.medaka.fasta; fi
  """
  stub:
  "touch ${meta.id}.medaka.fasta"
}
