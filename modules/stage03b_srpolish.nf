// Stage 03b — short-read (Illumina) polishing of the Medaka consensus with
// Polypolish. This is the HYBRID branch and becomes the default once short reads
// are present (--hybrid auto). Polypolish requires ALL alignments per read
// (bwa mem -a) and an insert-size filter, else it under-/mis-polishes. The mode
// reported here is HONEST: if Illumina is absent it is 'ont_only'; if the
// short-read pass is attempted but fails, it is 'hybrid_failed_ont_fallback'
// (never silently claimed as a clean hybrid). resistance_confidence rides
// downstream to flag ONT-only calls as provisional (homopolymer indels).
//
// Short-read-only isolates (assembly_mode == 'shortread') arrive here with a SPAdes
// assembly already built from the Illumina reads: there is nothing to hybrid-polish
// and no homopolymer-indel problem, so the assembly passes through unchanged with
// mode='illumina_only' and resistance_confidence='high'.
process SRPOLISH {
  tag { meta.id }
  label 'srpolish'
  publishDir { "${params.outdir}/${meta.id}/03_polish" }, mode: 'copy'
  input:  tuple val(meta), path(medaka), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.polished.fasta"), emit: assembly
          tuple val(meta), path("${meta.id}.polish.json"),    emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def is_shortread = meta.assembly_mode == 'shortread'
  def do_hybrid = !is_shortread && (params.hybrid == 'on' || (params.hybrid == 'auto' && has_illumina))
  if (is_shortread)
    """
    source "${projectDir}/bin/ff_status.sh"; ff_init polish "${meta.id}" ${meta.id}.polish.json
    ff_skip polypolish "Illumina-only assembly: short-read base accuracy needs no polishing"
    cp ${medaka} ${meta.id}.polished.fasta
    printf '{"sample":"%s","stage":"polish","mode":"illumina_only","short_read_polisher":"none","assembler":"%s","polypolish_bp_changed":"NA","resistance_confidence":"high"}\\n' \\
      "${meta.id}" "${params.sr_assembler}" > ${meta.id}.polish.json
    ff_finalize
    """
  else if (do_hybrid)
    """
    source "${projectDir}/bin/ff_status.sh"; ff_init polish "${meta.id}" ${meta.id}.polish.json
    ff_version polypolish -- polypolish --version
    ff_version bwa -- bwa
    # Each step is optional: a failure falls back to the medaka assembly, the mode records it
    # and the stage is marked partial (resistance confidence drops to provisional).
    MODE=hybrid_failed_ont_fallback; CONF=provisional_ont_only; CHANGES=NA; OK=1
    ff_run bwa_index --optional -- bash -c "bwa index ${medaka} 2> bwa_index.log"
    if [ "\$FF_RC" -ne 0 ]; then OK=0; fi
    if [ "\$OK" -eq 1 ]; then
      ff_run bwa_mem --optional -- bash -c "bwa mem -t ${task.cpus} -a ${medaka} ${r1} > r1.sam 2>> bwa.log && bwa mem -t ${task.cpus} -a ${medaka} ${r2} > r2.sam 2>> bwa.log"
      if [ "\$FF_RC" -ne 0 ]; then OK=0; fi
    fi
    if [ "\$OK" -eq 1 ]; then
      ff_run polypolish --optional -- bash -c "polypolish filter --in1 r1.sam --in2 r2.sam --out1 f1.sam --out2 f2.sam 2> polypolish_filter.log && polypolish polish ${medaka} f1.sam f2.sam > ${meta.id}.polished.fasta 2> polypolish.log && [ -s ${meta.id}.polished.fasta ]"
      if [ "\$FF_RC" -ne 0 ]; then OK=0; fi
    fi
    if [ "\$OK" -eq 1 ]; then
      MODE=hybrid; CONF=high
      CHANGES=\$(grep -oiE '[0-9,]+ bp changed|changed [0-9,]+' polypolish.log | grep -oE '[0-9,]+' | head -1 | tr -d ',' || echo NA)
    else
      cp ${medaka} ${meta.id}.polished.fasta
    fi
    printf '{"sample":"%s","stage":"polish","mode":"%s","short_read_polisher":"polypolish","polypolish_bp_changed":"%s","resistance_confidence":"%s"}\\n' \\
      "${meta.id}" "\$MODE" "\$CHANGES" "\$CONF" > ${meta.id}.polish.json
    ff_finalize
    """
  else
    """
    source "${projectDir}/bin/ff_status.sh"; ff_init polish "${meta.id}" ${meta.id}.polish.json
    ff_skip polypolish "no Illumina reads for this isolate (ONT-only)"
    cp ${medaka} ${meta.id}.polished.fasta
    printf '{"sample":"%s","stage":"polish","mode":"ont_only","short_read_polisher":"none","polypolish_bp_changed":"NA","resistance_confidence":"provisional_ont_only"}\\n' \\
      "${meta.id}" > ${meta.id}.polish.json
    ff_finalize
    """
  stub:
  "touch ${meta.id}.polished.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"polish\"}' > ${meta.id}.polish.json"
}
