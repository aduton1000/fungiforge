// Stage 02 — de-novo assembly (Flye --nano-hq default) + optional purge_dups to
// collapse heterozygous haplotigs. Emits a single nuclear+organelle draft.
process ASSEMBLE {
  tag { meta.id }
  label 'assembly'
  publishDir { "${params.outdir}/${meta.id}/02_assembly" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.assembly.fasta"), emit: assembly
  script:
  """
  case ${params.assembler} in
    flye)  flye --nano-hq ${ont} --out-dir flye --threads ${task.cpus} ; cp flye/assembly.fasta asm.fasta ;;
    raven) raven -t ${task.cpus} ${ont} > asm.fasta ;;
    canu)  canu -p ${meta.id} -d canu genomeSize=35m -nanopore ${ont} useGrid=false maxThreads=${task.cpus} ; cp canu/${meta.id}.contigs.fasta asm.fasta ;;
  esac
  if [ "${params.purge_dups}" = "true" ]; then
    minimap2 -t ${task.cpus} -xmap-ont asm.fasta ${ont} | gzip > aln.paf.gz || true
    ( purge_dups -h >/dev/null 2>&1 && \\
      pbcstat aln.paf.gz && calcuts PB.stat > cutoffs 2>calcults.log && \\
      split_fa asm.fasta > asm.split && minimap2 -t ${task.cpus} -xasm5 -DP asm.split asm.split | gzip > self.paf.gz && \\
      purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed 2>purge.log && \\
      get_seqs -e dups.bed asm.fasta && cp purged.fa ${meta.id}.assembly.fasta ) || cp asm.fasta ${meta.id}.assembly.fasta
  else
    cp asm.fasta ${meta.id}.assembly.fasta
  fi
  """
  stub:
  "touch ${meta.id}.assembly.fasta"
}
