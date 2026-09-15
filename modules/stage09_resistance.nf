// Stage 09 — antifungal resistance (bespoke). Extracts & aligns the target genes
// (CYP51A/ERG11, CYP51B, FKS1/2 HS1/HS2, FUR1, ERG3, HMG1, TAC1/MRR1/UPC2/PDR1,
// efflux CDR1/2/MDR1) against the curated FungAMR/MARDy allele+mutation panel, and
// runs the dedicated A. fumigatus cyp51A promoter TR34/TR46 detector. Calls carry
// a confidence flag (homopolymer / ONT-only-polish sensitive).
process RESISTANCE {
  tag { meta.id }
  label 'resistance'
  publishDir { "${params.outdir}/${meta.id}/09_resistance" }, mode: 'copy'
  input:  tuple val(meta), path(proteins), path(species), path(nuclear), path(gbk), path(polish_json)
  output: tuple val(meta), path("${meta.id}.resistance.json"), emit: json
  script:
  """
  source "${projectDir}/bin/ff_status.sh"; ff_init resistance "${meta.id}" ${meta.id}.resistance.json
  ff_version python -- python3 --version
  ff_version biopython -- python3 -c "import Bio; print(Bio.__version__)"
  ff_version hmmer -- bash -c "hmmsearch -h | grep -m1 HMMER"
  ff_version diamond -- diamond version
  ff_version minimap2 -- minimap2 --version
  # Carry the polishing mode through so ONT-only calls are flagged provisional while
  # hybrid-polished calls are high-confidence (the polish.json rides in from Stage 03b).
  MODE=\$(grep -o '"mode"[^,]*' ${polish_json} 2>/dev/null | head -1 | sed -E 's/.*: *"?([a-z_]+)"?.*/\\1/')
  case "\$MODE" in hybrid) PM=hybrid;; illumina_only) PM=illumina_only;; ont_only) PM=ont_only;; *) PM=unknown;; esac
  ff_run af_resistance -- python3 ${projectDir}/bin/af_resistance.py \\
      --sample "${meta.id}" --proteins ${proteins} --species ${species} \\
      --assembly ${nuclear} --gbk ${gbk} --polish-mode "\$PM" \\
      --panel "${params.af_panel}" --data-dir "${params.data_dir ?: ''}" \\
      --out ${meta.id}.resistance.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"resistance\"}' > ${meta.id}.resistance.json"
}
