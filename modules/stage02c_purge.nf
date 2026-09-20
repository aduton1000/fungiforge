// Stage 02c — haplotig purging (purge_dups, W4.1/L30). Long-read assemblies of heterozygous
// isolates carry duplicated haplotigs; purge_dups collapses them using the read depth of the
// ONT reads mapped back to the draft. It runs in the base image because the assembler images
// (staphb/flye, staphb/spades) carry neither purge_dups nor minimap2 — running it inside
// stage 02 exited 127 on every isolate and made the assembly stage permanently `partial`.
// Short-read assemblies skip it (SPAdes already resolves bubbles). --purge_dups false disables it.
// Organelle contigs are restored afterwards: their read depth is far above the nuclear peak, so
// purge_dups drops them as collapsed duplication (it removed the CEA10 mitogenome outright).
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
      # purge_dups judges by read depth against the nuclear peak, and a mitochondrial genome sits
      # far above it, so it is dropped as collapsed duplication. On CEA10 that removed the 30.7 kb,
      # 25.5 % GC mitogenome and emptied every mitochondrial column downstream. Put back only the
      # contigs that carry core mitochondrial genes; genuine haplotigs stay purged.
      if [ -s "${params.mito_proteins}" ]; then
        ff_run makeblastdb_pre --optional -- bash -c "makeblastdb -in ${assembly} -dbtype nucl -out preasm_db > makeblastdb.log 2>&1"
        if [ "\$FF_RC" -eq 0 ]; then
          ff_run tblastn_mito_pre --optional -- bash -c "tblastn -query '${params.mito_proteins}' -db preasm_db -evalue 1e-10 -max_target_seqs 50 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > mito_hits.tsv 2>tblastn.log"
        fi
        if [ "\$FF_RC" -eq 0 ]; then
          ff_run protect_organelle --optional -- protect_organelle.py --assembly ${assembly} \\
              --purged ${meta.id}.purged.fasta --tblastn mito_hits.tsv --max-len ${params.mito_max_len} \\
              --out restored.fasta --json protect_organelle.json
          if [ "\$FF_RC" -eq 0 ] && [ -s restored.fasta ]; then mv restored.fasta ${meta.id}.purged.fasta; fi
        fi
      else
        ff_skip protect_organelle "no --mito_proteins: cannot tell an organelle contig from a haplotig, nothing restored"
      fi
      N_PURGED=\$(( \$(grep -c '^>' ${assembly}) - \$(grep -c '^>' ${meta.id}.purged.fasta) ))
    fi
  fi
  [ -s ${meta.id}.purged.fasta ] || cp ${assembly} ${meta.id}.purged.fasta
  N_RESTORED=\$(python3 -c "import json,sys;print(json.load(open('protect_organelle.json')).get('n_restored',0))" 2>/dev/null || echo 0)
  printf '{"sample":"%s","stage":"purge","purge_dups_applied":%s,"n_contigs_before":%s,"n_contigs_after":%s,"n_haplotigs_removed":%s,"n_organelle_restored":%s}\\n' \\
    "${meta.id}" "\$PURGED" "\$(grep -c '^>' ${assembly})" "\$(grep -c '^>' ${meta.id}.purged.fasta)" "\$N_PURGED" "\${N_RESTORED:-0}" > ${meta.id}.purge.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.purged.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"purge\"}' > ${meta.id}.purge.json"
}
