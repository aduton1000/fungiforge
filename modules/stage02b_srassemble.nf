// Stage 02b — short-read (Illumina-only) de-novo assembly. Used for isolates with
// no ONT reads (meta.assembly_mode == 'shortread'). SPAdes --isolate by default
// (fungal-genome friendly, built-in read error-correction), MEGAHIT as a lean
// alternative. Emits the same ${meta.id}.assembly.fasta contract as Stage 02 so the
// two assembly paths merge transparently in the subworkflow. No purge_dups here
// (that is a long-read haplotig step); short-read contigs are already collapsed.
// Illumina base accuracy means there is no homopolymer-indel problem, so these
// assemblies carry HIGH-confidence resistance calls (see Stage 03b / af_resistance).
process SR_ASSEMBLE {
  tag { meta.id }
  label 'srassembly'
  publishDir { "${params.outdir}/${meta.id}/02_assembly" }, mode: 'copy'
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.assembly.fasta"), emit: assembly
  script:
  def mem_gb = (task.memory ? task.memory.toGiga() : 64)
  """
  case ${params.sr_assembler} in
    spades)
      spades.py --isolate -1 ${r1} -2 ${r2} -o spades -t ${task.cpus} -m ${mem_gb}
      cp spades/scaffolds.fasta ${meta.id}.assembly.fasta 2>/dev/null \\
        || cp spades/contigs.fasta ${meta.id}.assembly.fasta
      ;;
    megahit)
      megahit -1 ${r1} -2 ${r2} -o megahit -t ${task.cpus}
      cp megahit/final.contigs.fa ${meta.id}.assembly.fasta
      ;;
    *)
      echo "Unknown --sr_assembler '${params.sr_assembler}' (use spades|megahit)" >&2; exit 1 ;;
  esac
  """
  stub:
  "touch ${meta.id}.assembly.fasta"
}
