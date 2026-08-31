// Stage 00 — optional Dorado basecalling (pod5 dir -> ONT FASTQ). Runs native
// arm64 Dorado on the host (fast); skipped unless --basecall.
process BASECALL {
  tag { meta.id }
  label 'basecall'
  publishDir { "${params.outdir}/${meta.id}/00_basecall" }, mode: 'copy'
  input:  tuple val(meta), path(pod5), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.ont.fastq.gz"), path(r1), path(r2), emit: reads
  script:
  """
  dorado basecaller ${params.dorado_duplex ? 'duplex' : params.dorado_model} ${pod5} \\
      --emit-fastq ${params.dorado_duplex ? '' : ''} > ${meta.id}.ont.fastq
  gzip -f ${meta.id}.ont.fastq
  """
  stub:
  "touch ${meta.id}.ont.fastq.gz"
}
