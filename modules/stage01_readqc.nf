// Stage 01 — read QC + filtering (mode-aware). Long-read isolates: NanoPlot
// (raw+filtered) + chopper on ONT, plus fastp on Illumina when present (hybrid).
// Short-read-only isolates (assembly_mode == 'shortread'): fastp on Illumina only,
// and an empty ONT placeholder is emitted to keep the downstream tuple contract.
// Emits the FILTERED ONT reads and the TRIMMED Illumina reads (W2.2: assembly and polishing
// now use fastp's output, not the raw pairs) plus a QC JSON with read statistics
// (fastp / NanoPlot numbers via bin/read_stats.py).
process READ_QC {
  tag { meta.id }
  label 'readqc'
  publishDir { "${params.outdir}/${meta.id}/01_readqc" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.ont.filt.fastq.gz"),
                          path("${r1.name == 'NO_R1' ? 'NO_R1' : meta.id + '_R1.trim.fastq.gz'}"),
                          path("${r2.name == 'NO_R2' ? 'NO_R2' : meta.id + '_R2.trim.fastq.gz'}"), emit: reads
          tuple val(meta), path("${meta.id}.readqc.json"),                                          emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def is_shortread = meta.assembly_mode == 'shortread'
  def platform     = is_shortread ? 'illumina' : (has_illumina ? 'hybrid' : 'ont')
  def fastp        = "ff_run fastp -- fastp -i ${r1} -I ${r2} -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz --json fastp.json --html fastp.html --thread ${task.cpus}"
  if (is_shortread)
    """
    source ff_status.sh; ff_init readqc "${meta.id}" ${meta.id}.readqc.json
    ff_version fastp -- fastp --version
    ${fastp}
    : | gzip > ${meta.id}.ont.filt.fastq.gz   # no ONT for this isolate — empty placeholder
    printf '{"sample":"%s","stage":"readqc","platform":"%s","ont_filtered":null,"illumina":true}\\n' \\
      "${meta.id}" "${platform}" > ${meta.id}.readqc.json
    ff_run read_stats -- read_stats.py --sample "${meta.id}" --platform ${platform} --fastp fastp.json \\
        --merge-into ${meta.id}.readqc.json --out ${meta.id}.readqc.json
    ff_finalize
    """
  else
    """
    source ff_status.sh; ff_init readqc "${meta.id}" ${meta.id}.readqc.json
    ff_version chopper -- chopper --version
    ff_version nanoplot -- NanoPlot --version
    ff_version fastp -- fastp --version
    # --plots dot --no_static: NanoPlot 1.47's 2-D density plot needs a plotly API removed in plotly 6, and
    # static PNG export needs a Chrome-based Kaleido; neither is in the image, and only NanoStats.txt is consumed.
    ff_run nanoplot_raw --optional -- NanoPlot --fastq ${ont} -o nanoplot_raw --prefix raw --N50 --plots dot --no_static
    ff_run chopper -- bash -o pipefail -c "chopper -q ${params.ont_min_qual} -l ${params.ont_min_len} -i ${ont} 2> chopper.log | gzip > ${meta.id}.ont.filt.fastq.gz"
    ff_run nanoplot_filt --optional -- NanoPlot --fastq ${meta.id}.ont.filt.fastq.gz -o nanoplot_filt --prefix filt --N50 --plots dot --no_static
    ${ has_illumina ? fastp : "ff_skip fastp 'no Illumina reads for this isolate'" }
    printf '{"sample":"%s","stage":"readqc","platform":"%s","ont_filtered":"%s","illumina":%s}\\n' \\
      "${meta.id}" "${platform}" "${meta.id}.ont.filt.fastq.gz" "${has_illumina.toString()}" > ${meta.id}.readqc.json
    ff_run read_stats -- read_stats.py --sample "${meta.id}" --platform ${platform} ${ has_illumina ? '--fastp fastp.json' : '' } \\
        --nanostats-raw nanoplot_raw/rawNanoStats.txt --nanostats-filt nanoplot_filt/filtNanoStats.txt \\
        --merge-into ${meta.id}.readqc.json --out ${meta.id}.readqc.json
    ff_finalize
    """
  stub:
  def r1o = r1.name == 'NO_R1' ? '' : "touch ${meta.id}_R1.trim.fastq.gz ${meta.id}_R2.trim.fastq.gz;"
  "touch ${meta.id}.ont.filt.fastq.gz; ${r1o} echo '{\"sample\":\"${meta.id}\",\"stage\":\"readqc\"}' > ${meta.id}.readqc.json"
}
