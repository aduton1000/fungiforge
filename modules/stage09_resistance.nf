// Stage 09 — antifungal resistance (bespoke). Extracts & aligns the target genes
// (CYP51A/ERG11, CYP51B, FKS1/2 HS1/HS2, FUR1, ERG3, HMG1, TAC1/MRR1/UPC2/PDR1,
// efflux CDR1/2/MDR1) against the curated FungAMR/MARDy allele+mutation panel, and
// runs the dedicated A. fumigatus cyp51A promoter TR34/TR46 detector. Calls carry
// a confidence flag (homopolymer / ONT-only-polish sensitive).
process RESISTANCE {
  tag { meta.id }
  label 'resistance'
  publishDir { "${params.outdir}/${meta.id}/09_resistance" }, mode: 'copy'
  input:  tuple val(meta), path(proteins), path(species), path(nuclear), path(gbk)
  output: tuple val(meta), path("${meta.id}.resistance.json"), emit: json
  script:
  """
  python3 ${projectDir}/bin/af_resistance.py \\
      --sample "${meta.id}" --proteins ${proteins} --species ${species} \\
      --assembly ${nuclear} --gbk ${gbk} \\
      --panel "${params.af_panel}" --data-dir "${params.data_dir ?: ''}" \\
      --out ${meta.id}.resistance.json || echo '{"sample":"${meta.id}","stage":"resistance","calls":[]}' > ${meta.id}.resistance.json
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"resistance\"}' > ${meta.id}.resistance.json"
}
