// Stage 02c — haplotig purging (purge_dups, W4.1/L30). Long-read assemblies of heterozygous
// isolates carry duplicated haplotigs; purge_dups collapses them using the read depth of the
// ONT reads mapped back to the draft. It runs in the base image because the assembler images
// (staphb/flye, staphb/spades) carry neither purge_dups nor minimap2 — running it inside
// stage 02 exited 127 on every isolate and made the assembly stage permanently `partial`.
// Short-read assemblies skip it (SPAdes already resolves bubbles). --purge_dups false disables it.
process PURGE_DUPS {
  tag { meta.id }
  label 'decontam'
  publishDir { "${params.outdir}/${meta.id}/02_assembly" }, mode: 'copy', pattern: "*.{json,bed}"
  input:  tuple val(meta), path(assembly), path(ont)
  output: tuple val(meta), path("${meta.id}.purged.fasta"), emit: assembly
          tuple val(meta), path("${meta.id}.purge.json"),   emit: json
  script:
  """
  source ff_status.sh; ff_init purge "${meta.id}" ${meta.id}.purge.json --best-effort
  ff_version purge_dups -- bash -c "purge_dups 2>&1 | grep -m1 -i version || echo unknown"
  ff_version minimap2 -- minimap2 --version
  PURGED=false; N_PURGED=0
  if [ "${params.purge_dups}" != "true" ]; then
    ff_skip purge_dups "disabled (--purge_dups false)"
  elif [ ! -s ${ont} ]; then
    ff_skip purge_dups "no ONT reads for this isolate"
  else
    ff_run purge_dups --optional -- bash -o pipefail -c "minimap2 -t ${task.cpus} -x map-ont ${assembly} ${ont} 2>minimap2.log | gzip > aln.paf.gz \\
         && pbcstat aln.paf.gz && calcuts PB.stat > cutoffs 2>calcuts.log \\
         && split_fa ${assembly} > asm.split && minimap2 -t ${task.cpus} -x asm5 -DP asm.split asm.split 2>>minimap2.log | gzip > self.paf.gz \\
         && purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed 2>purge.log \\
         && get_seqs -e dups.bed ${assembly} && [ -s purged.fa ]"
    if [ "\$FF_RC" -eq 0 ]; then
      cp purged.fa ${meta.id}.purged.fasta; PURGED=true
      cp dups.bed ${meta.id}.dups.bed 2>/dev/null || true
      N_PURGED=\$(( \$(grep -c '^>' ${assembly}) - \$(grep -c '^>' ${meta.id}.purged.fasta) ))
    fi
  fi
  [ -s ${meta.id}.purged.fasta ] || cp ${assembly} ${meta.id}.purged.fasta
  printf '{"sample":"%s","stage":"purge","purge_dups_applied":%s,"n_contigs_before":%s,"n_contigs_after":%s,"n_haplotigs_removed":%s}\\n' \\
    "${meta.id}" "\$PURGED" "\$(grep -c '^>' ${assembly})" "\$(grep -c '^>' ${meta.id}.purged.fasta)" "\$N_PURGED" > ${meta.id}.purge.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.purged.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"purge\"}' > ${meta.id}.purge.json"
}
