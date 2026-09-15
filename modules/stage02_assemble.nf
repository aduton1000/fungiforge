// Stage 02 — de-novo assembly (Flye --nano-hq default) + optional purge_dups to
// collapse heterozygous haplotigs. Emits a single nuclear+organelle draft.
process ASSEMBLE {
  tag { meta.id }
  label 'assembly'
  publishDir { "${params.outdir}/${meta.id}/02_assembly" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.assembly.fasta"), emit: assembly
          tuple val(meta), path("${meta.id}.assemble.json"),   emit: json
  script:
  """
  source ff_status.sh; ff_init assemble "${meta.id}" ${meta.id}.assemble.json
  ff_version ${params.assembler} -- ${params.assembler == "flye" ? "flye --version" : params.assembler == "raven" ? "raven --version" : "canu -version"}
  ff_version minimap2 -- minimap2 --version
  ff_version purge_dups -- purge_dups -h
  case ${params.assembler} in
    flye)  ff_run flye -- flye --nano-hq ${ont} --out-dir flye --threads ${task.cpus}; cp flye/assembly.fasta asm.fasta ;;
    raven) ff_run raven -- bash -c "raven -t ${task.cpus} ${ont} > asm.fasta" ;;
    canu)  ff_run canu -- canu -p ${meta.id} -d canu genomeSize=35m -nanopore ${ont} useGrid=false maxThreads=${task.cpus}; cp canu/${meta.id}.contigs.fasta asm.fasta ;;
    *)     echo "Unknown --assembler '${params.assembler}' (use flye|canu|raven)" >&2; exit 1 ;;
  esac
  PURGED=false
  if [ "${params.purge_dups}" = "true" ]; then
    # optional: a purge_dups failure falls back to the unpurged assembly and marks the stage partial
    ff_run purge_dups --optional -- bash -o pipefail -c "minimap2 -t ${task.cpus} -xmap-ont asm.fasta ${ont} | gzip > aln.paf.gz \\
         && pbcstat aln.paf.gz && calcuts PB.stat > cutoffs 2>calcults.log \\
         && split_fa asm.fasta > asm.split && minimap2 -t ${task.cpus} -xasm5 -DP asm.split asm.split | gzip > self.paf.gz \\
         && purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed 2>purge.log \\
         && get_seqs -e dups.bed asm.fasta && [ -s purged.fa ]"
    if [ "\$FF_RC" -eq 0 ]; then
      cp purged.fa ${meta.id}.assembly.fasta; PURGED=true
    else
      cp asm.fasta ${meta.id}.assembly.fasta
    fi
  else
    ff_skip purge_dups "disabled (--purge_dups false)"
    cp asm.fasta ${meta.id}.assembly.fasta
  fi
  printf '{"sample":"%s","stage":"assemble","assembler":"%s","purge_dups_applied":%s,"n_contigs":%s}\\n' \\
    "${meta.id}" "${params.assembler}" "\$PURGED" "\$(grep -c '^>' ${meta.id}.assembly.fasta)" > ${meta.id}.assemble.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.assembly.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"assemble\"}' > ${meta.id}.assemble.json"
}
