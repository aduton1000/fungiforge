// Stage 01 — read QC + filtering (mode-aware). Long-read isolates: NanoPlot
// (raw+filtered) + chopper on ONT, plus fastp on Illumina when present (hybrid).
// Short-read-only isolates (assembly_mode == 'shortread'): fastp on Illumina only,
// and an empty ONT placeholder is emitted to keep the downstream tuple contract.
// Emits filtered/passed reads + a QC result.json carrying the platform.
process READ_QC {
  tag { meta.id }
  label 'readqc'
  publishDir { "${params.outdir}/${meta.id}/01_readqc" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.ont.filt.fastq.gz"), path(r1), path(r2), emit: reads
          tuple val(meta), path("${meta.id}.readqc.json"),                           emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def is_shortread = meta.assembly_mode == 'shortread'
  def platform     = is_shortread ? 'illumina' : (has_illumina ? 'hybrid' : 'ont')
  if (is_shortread)
    """
    fastp -i ${r1} -I ${r2} -o r1.fp.fq.gz -O r2.fp.fq.gz --json fastp.json --thread ${task.cpus} || true
    : | gzip > ${meta.id}.ont.filt.fastq.gz   # no ONT for this isolate — empty placeholder
    printf '{"sample":"%s","stage":"readqc","platform":"%s","ont_filtered":null,"illumina":true}\\n' \\
      "${meta.id}" "${platform}" > ${meta.id}.readqc.json
    """
  else
    """
    NanoPlot --fastq ${ont} -o nanoplot_raw --prefix raw    --N50 || true
    chopper -q ${params.ont_min_qual} -l ${params.ont_min_len} -i ${ont} 2> chopper.log \\
        | gzip > ${meta.id}.ont.filt.fastq.gz
    NanoPlot --fastq ${meta.id}.ont.filt.fastq.gz -o nanoplot_filt --prefix filt --N50 || true
    ${ has_illumina ? "fastp -i ${r1} -I ${r2} -o r1.fp.fq.gz -O r2.fp.fq.gz --json fastp.json --thread ${task.cpus} || true" : "" }
    printf '{"sample":"%s","stage":"readqc","platform":"%s","ont_filtered":"%s","illumina":%s}\\n' \\
      "${meta.id}" "${platform}" "${meta.id}.ont.filt.fastq.gz" "${has_illumina.toString()}" > ${meta.id}.readqc.json
    """
  stub:
  "touch ${meta.id}.ont.filt.fastq.gz ${meta.id}.readqc.json"
}
