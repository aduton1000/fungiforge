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
          tuple val(meta), path("${meta.id}.assemble.json"),   emit: json
  script:
  def mem_gb = (task.memory ? task.memory.toGiga() : 64)
  """
  source "${projectDir}/bin/ff_status.sh"; ff_init assemble "${meta.id}" ${meta.id}.assemble.json
  case ${params.sr_assembler} in
    spades)
      ff_run spades -- spades.py --isolate -1 ${r1} -2 ${r2} -o spades -t ${task.cpus} -m ${mem_gb}
      if [ -s spades/scaffolds.fasta ]; then cp spades/scaffolds.fasta raw.fasta; else cp spades/contigs.fasta raw.fasta; fi
      ;;
    megahit)
      ff_run megahit -- megahit -1 ${r1} -2 ${r2} -o megahit -t ${task.cpus}
      cp megahit/final.contigs.fa raw.fasta
      ;;
    *)
      echo "Unknown --sr_assembler '${params.sr_assembler}' (use spades|megahit)" >&2; exit 1 ;;
  esac
  # Keep contigs >= --min_contig_len with >= 4 distinct bases. Short-read assemblers emit
  # thousands of tiny fragments and the odd homopolymer stub (e.g. 78 bp of G); funannotate
  # refuses an assembly containing the latter, and the former only inflate contig counts.
  ff_run contig_filter -- awk -v min=${params.min_contig_len} 'BEGIN{RS=">"; ORS=""}
    NR>1 { n=index(\$0,"\\n"); hdr=substr(\$0,1,n-1); seq=substr(\$0,n+1); gsub(/\\n/,"",seq)
           if (length(seq) < min) next
           s=seq; a=gsub(/[Aa]/,"",s); s=seq; c=gsub(/[Cc]/,"",s); s=seq; g=gsub(/[Gg]/,"",s); s=seq; t=gsub(/[Tt]/,"",s)
           if ((a>0)+(c>0)+(g>0)+(t>0) < 4) next
           printf(">%s\\n%s\\n", hdr, seq) }' raw.fasta > ${meta.id}.assembly.fasta
  KEPT=\$(grep -c '^>' ${meta.id}.assembly.fasta); RAW=\$(grep -c '^>' raw.fasta)
  echo "[sr_assemble] kept \$KEPT of \$RAW contigs (>= ${params.min_contig_len} bp, >= 4 bases)"
  printf '{"sample":"%s","stage":"assemble","assembler":"%s","n_contigs_raw":%s,"n_contigs":%s,"min_contig_len":%s}\\n' \\
    "${meta.id}" "${params.sr_assembler}" "\$RAW" "\$KEPT" "${params.min_contig_len}" > ${meta.id}.assemble.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.assembly.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"assemble\"}' > ${meta.id}.assemble.json"
}
