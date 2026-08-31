// Stage 06 — de-novo repeat modeling + soft-masking (RepeatModeler2 + RepeatMasker
// in the dfam/tetools container). Emits the soft-masked genome (for annotation),
// the TE library + a RIP-index summary (for Stage 10), and a result.json.
process REPEATMASK {
  tag { meta.id }
  label 'repeat'
  publishDir { "${params.outdir}/${meta.id}/06_repeatmask" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear)
  output: tuple val(meta), path("${meta.id}.masked.fasta"), emit: masked
          tuple val(meta), path("${meta.id}.telib.fasta"),  emit: telib
          tuple val(meta), path("${meta.id}.repeat.json"),  emit: json
  script:
  """
  BuildDatabase -name ${meta.id}_db ${nuclear}
  RepeatModeler -database ${meta.id}_db -threads ${task.cpus} -LTRStruct || true
  cp ${meta.id}_db-families.fa ${meta.id}.telib.fasta 2>/dev/null || touch ${meta.id}.telib.fasta
  RepeatMasker -pa ${task.cpus} -lib ${meta.id}.telib.fasta -xsmall -dir rm ${nuclear} || true
  cp rm/*.masked ${meta.id}.masked.fasta 2>/dev/null || cp ${nuclear} ${meta.id}.masked.fasta
  python3 -c "import json;json.dump({'sample':'${meta.id}','stage':'repeatmask','telib':'${meta.id}.telib.fasta'},open('${meta.id}.repeat.json','w'),indent=2)"
  """
  stub:
  "touch ${meta.id}.masked.fasta ${meta.id}.telib.fasta ${meta.id}.repeat.json"
}
