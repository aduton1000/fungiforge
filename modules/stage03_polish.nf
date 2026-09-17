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
  # Model choice: Dorado writes basecall_model_version_id=... into every read header and
  # medaka resolves its model from that by itself; reads without it (other basecallers, public
  # data) use --medaka_model (default: the --dorado_model string), which must be a model this
  # medaka knows — an unknown model is a failed stage, not a silent default.
  HDR_MODEL=\$( (zcat -f ${ont} 2>/dev/null || cat ${ont}) | head -1 | grep -oE 'basecall_model_version_id=[^[:space:]]+' | cut -d= -f2 || true)
  if [ -n "\$HDR_MODEL" ]; then
    MODEL="auto(\$HDR_MODEL)"; MODEL_ARG=""
  else
    MODEL="${params.medaka_model ?: params.dorado_model}"; MODEL_ARG="-m \$MODEL"
    ff_run medaka_model_check -- bash -c "medaka tools list_models | tr ' ,' '\\n\\n' | grep -qx '\$MODEL'"
  fi
  ff_run medaka -- medaka_consensus -i ${ont} -d ${assembly} -o medaka -t ${task.cpus} \$MODEL_ARG
  ff_run medaka_output -- test -s medaka/consensus.fasta
  cp medaka/consensus.fasta ${meta.id}.medaka.fasta
  printf '{"sample":"%s","stage":"medaka","model":"%s"}\\n' "${meta.id}" "\$MODEL" > ${meta.id}.medaka.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.medaka.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"medaka\"}' > ${meta.id}.medaka.json"
}
