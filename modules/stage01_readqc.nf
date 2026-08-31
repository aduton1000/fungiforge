// Stage 01 — read QC + filtering. ONT: NanoPlot (raw+filtered) + chopper.
// Illumina (if present): fastp. Emits filtered reads + a QC result.json.
process READ_QC {
  tag { meta.id }
  label 'readqc'
  publishDir { "${params.outdir}/${meta.id}/01_readqc" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.ont.filt.fastq.gz"), path(r1), path(r2), emit: reads
          tuple val(meta), path("${meta.id}.readqc.json"),                           emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  """
  NanoPlot --fastq ${ont} -o nanoplot_raw --prefix raw    --N50 || true
  chopper -q ${params.ont_min_qual} -l ${params.ont_min_len} -i ${ont} 2> chopper.log \\
      | gzip > ${meta.id}.ont.filt.fastq.gz
  NanoPlot --fastq ${meta.id}.ont.filt.fastq.gz -o nanoplot_filt --prefix filt --N50 || true
  ${ has_illumina ? "fastp -i ${r1} -I ${r2} -o r1.fp.fq.gz -O r2.fp.fq.gz --json fastp.json --thread ${task.cpus} || true" : "" }
  printf '{"sample":"%s","stage":"readqc","ont_filtered":"%s","illumina":%s}\\n' \\
    "${meta.id}" "${meta.id}.ont.filt.fastq.gz" "${has_illumina.toString()}" > ${meta.id}.readqc.json
  """
  stub:
  "touch ${meta.id}.ont.filt.fastq.gz ${meta.id}.readqc.json"
}
