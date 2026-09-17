// Stage 02 — de-novo assembly (Flye --nano-hq default). Emits a single nuclear+organelle draft.
// Haplotig purging is stage 02c (its own process in the base image): the assembler images carry
// neither purge_dups nor minimap2, so running it here exited 127 on every long-read isolate and
// left the stage `partial` (L30).
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
  case ${params.assembler} in
    flye)  ff_run flye -- flye --nano-hq ${ont} --out-dir flye --threads ${task.cpus}; cp flye/assembly.fasta asm.fasta ;;
    raven) ff_run raven -- bash -c "raven -t ${task.cpus} ${ont} > asm.fasta" ;;
    canu)  ff_run canu -- canu -p ${meta.id} -d canu genomeSize=35m -nanopore ${ont} useGrid=false maxThreads=${task.cpus}; cp canu/${meta.id}.contigs.fasta asm.fasta ;;
    *)     echo "Unknown --assembler '${params.assembler}' (use flye|canu|raven)" >&2; exit 1 ;;
  esac
  cp asm.fasta ${meta.id}.assembly.fasta
  printf '{"sample":"%s","stage":"assemble","assembler":"%s","purge_dups_applied":%s,"n_contigs":%s}\\n' \\
    "${meta.id}" "${params.assembler}" "\$PURGED" "\$(grep -c '^>' ${meta.id}.assembly.fasta)" > ${meta.id}.assemble.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.assembly.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"assemble\"}' > ${meta.id}.assemble.json"
}
