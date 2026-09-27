// Stage 00 — optional Dorado basecalling (pod5 dir -> ONT FASTQ). Runs native
// Dorado as a native host binary (fast); skipped unless --basecall.
process BASECALL {
  tag { meta.id }
  label 'basecall'
  publishDir { "${params.outdir}/${meta.id}/00_basecall" }, mode: 'copy'
  input:  tuple val(meta), path(pod5), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.ont.fastq.gz"), path(r1), path(r2), emit: reads
          tuple val(meta), path("${meta.id}.basecall.json"),  emit: json
  script:
  """
  source ff_status.sh; ff_init basecall "${meta.id}" ${meta.id}.basecall.json
  ff_version dorado -- dorado --version
  ff_run dorado -- bash -o pipefail -c "dorado basecaller ${params.dorado_duplex ? 'duplex' : params.dorado_model} ${pod5} --emit-fastq > ${meta.id}.ont.fastq"
  ff_run gzip -- gzip -f ${meta.id}.ont.fastq
  printf '{"sample":"%s","stage":"basecall","model":"%s","duplex":%s}\\n' "${meta.id}" "${params.dorado_model}" "${params.dorado_duplex}" > ${meta.id}.basecall.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.ont.fastq.gz; echo '{\"sample\":\"${meta.id}\",\"stage\":\"basecall\"}' > ${meta.id}.basecall.json"
}
