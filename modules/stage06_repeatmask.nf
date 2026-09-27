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
          tuple val(meta), path("${meta.id}.repeatmasker.tbl"), emit: tbl    // W2.8: RepeatMasker summary for te_percent (empty when not run)
  script:
  """
  source ff_status.sh; ff_init repeatmask "${meta.id}" ${meta.id}.repeat.json
  ff_version RepeatModeler -- RepeatModeler -version
  ff_version RepeatMasker -- RepeatMasker -v
  ff_run BuildDatabase -- BuildDatabase -name ${meta.id}_db ${nuclear}
  ff_run RepeatModeler -- RepeatModeler -database ${meta.id}_db -threads ${task.cpus} -LTRStruct
  MASKED=false; NFAM=0
  if [ -s ${meta.id}_db-families.fa ]; then
    cp ${meta.id}_db-families.fa ${meta.id}.telib.fasta
    NFAM=\$(grep -c '^>' ${meta.id}.telib.fasta)
    ff_run RepeatMasker -- RepeatMasker -pa ${task.cpus} -lib ${meta.id}.telib.fasta -xsmall -dir rm ${nuclear}
    if ls rm/*.masked >/dev/null 2>&1; then cp rm/*.masked ${meta.id}.masked.fasta; MASKED=true; else cp ${nuclear} ${meta.id}.masked.fasta; fi
    if ls rm/*.tbl >/dev/null 2>&1; then cp rm/*.tbl ${meta.id}.repeatmasker.tbl; fi
  else
    # RepeatModeler found no repeat families (possible for a small, repeat-poor genome):
    # nothing to mask — recorded, not hidden
    : > ${meta.id}.telib.fasta
    ff_skip RepeatMasker "RepeatModeler produced no repeat families; assembly left unmasked"
    cp ${nuclear} ${meta.id}.masked.fasta
  fi
  [ -f ${meta.id}.repeatmasker.tbl ] || : > ${meta.id}.repeatmasker.tbl
  printf '{"sample":"%s","stage":"repeatmask","telib":"%s.telib.fasta","n_families":%s,"masked":%s}\\n' \\
    "${meta.id}" "${meta.id}" "\$NFAM" "\$MASKED" > ${meta.id}.repeat.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.masked.fasta ${meta.id}.telib.fasta ${meta.id}.repeatmasker.tbl; echo '{\"sample\":\"${meta.id}\",\"stage\":\"repeatmask\"}' > ${meta.id}.repeat.json"
}
