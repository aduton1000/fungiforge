// Stage 09 — antifungal resistance (bespoke). Extracts & aligns the target genes
// (CYP51A/ERG11, CYP51B, FKS1/2 HS1/HS2, FUR1, ERG3, HMG1, TAC1/MRR1/UPC2/PDR1,
// efflux CDR1/2/MDR1) against the curated panel merged with the FungAMR catalogue (evidence
// tiers), runs the A. fumigatus cyp51A promoter TR34/TR46 detector, and (W2.4) re-genotypes
// every hotspot from the reads: a subsample of the QC'd reads is mapped to the GenBank records
// (the coordinate system of the annotation), giving per-hotspot allele frequency and zygosity,
// the TR site as a read-level insertion/deletion, and locus/genome depth as a copy-number
// signal. Calls carry a confidence combining polish mode, read support and evidence tier.
process RESISTANCE {
  tag { meta.id }
  label 'resistance'
  publishDir { "${params.outdir}/${meta.id}/09_resistance" }, mode: 'copy', pattern: "*.{json,bam,bai}"
  input:  tuple val(meta), path(proteins), path(species), path(nuclear), path(gbk), path(polish_json), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.resistance.json"), emit: json
          tuple val(meta), path("${meta.id}.reads.bam"), path("${meta.id}.reads.bam.bai"), emit: bam   // W2.6: reused by stage 13 (nQuire)
  script:
  def has_illumina = r1.name != 'NO_R1'
  def n_ilmn = params.genotype_reads
  def n_ont  = params.genotype_reads_ont
  def sub_cmd = has_illumina ? "seqkit head -n ${n_ilmn} ${r1} > sub_1.fq && seqkit head -n ${n_ilmn} ${r2} > sub_2.fq"
                             : "seqkit head -n ${n_ont} ${ont} > sub.fq"
  def map_cmd = has_illumina ? "minimap2 -ax sr -t ${task.cpus} --secondary=no records.fa sub_1.fq sub_2.fq"
                             : "minimap2 -ax map-ont -t ${task.cpus} --secondary=no records.fa sub.fq"
  """
  source ff_status.sh; ff_init resistance "${meta.id}" ${meta.id}.resistance.json
  ff_version python -- python3 --version
  ff_version biopython -- python3 -c "import Bio; print(Bio.__version__)"
  ff_version minimap2 -- minimap2 --version
  ff_version samtools -- samtools --version
  ff_version seqkit -- seqkit version
  # Carry the polishing mode through so ONT-only calls are flagged provisional while
  # hybrid-polished calls are high-confidence (the polish.json rides in from Stage 03b).
  MODE=\$(grep -o '"mode"[^,]*' ${polish_json} 2>/dev/null | head -1 | sed -E 's/.*: *"?([a-z_]+)"?.*/\\1/')
  case "\$MODE" in hybrid) PM=hybrid;; illumina_only) PM=illumina_only;; ont_only) PM=ont_only;; *) PM=unknown;; esac
  BAM_ARG=""
  if [ "${params.read_genotype}" = "true" ]; then
    ff_run gbk_to_fasta --optional -- gbk_to_fasta.py --gbk ${gbk} --out records.fa
    if [ "\$FF_RC" -eq 0 ] && [ -s records.fa ]; then
      ff_run subsample --optional -- bash -o pipefail -c "${sub_cmd}"
      ff_run map_reads --optional -- bash -o pipefail -c "${map_cmd} 2>minimap2.log | samtools sort -@ ${task.cpus} -m 1G -o ${meta.id}.reads.bam - && samtools index ${meta.id}.reads.bam"
      [ "\$FF_RC" -eq 0 ] && BAM_ARG="--bam ${meta.id}.reads.bam"
    fi
  else
    ff_skip map_reads "read-level genotyping disabled (--read_genotype false)"
  fi
  # stage 13 joins on the BAM: an empty placeholder pair keeps the contract when mapping did not run
  [ -s ${meta.id}.reads.bam ] || { : > ${meta.id}.reads.bam; : > ${meta.id}.reads.bam.bai; }
  ff_run af_resistance -- af_resistance.py \\
      --sample "${meta.id}" --proteins ${proteins} --species ${species} \\
      --assembly ${nuclear} --gbk ${gbk} --polish-mode "\$PM" \$BAM_ARG \\
      --min-reads ${params.genotype_min_reads} --min-frac ${params.genotype_min_frac} \\
      --panel "${params.af_panel}" --data-dir "${params.data_dir ?: ''}" \\
      --out ${meta.id}.resistance.json
  ff_finalize
  """
  stub:
  ": > ${meta.id}.reads.bam; : > ${meta.id}.reads.bam.bai; echo '{\"sample\":\"${meta.id}\",\"stage\":\"resistance\"}' > ${meta.id}.resistance.json"
}
