// Stage 13 — eukaryote-specific extras. Ploidy/heterozygosity (nQuire + Smudgeplot
// + GenomeScope2), mating-type locus (MAT1-1/MAT1-2 idiomorphs), secretome
// (SignalP6 + DeepTMHMM + EffectorP3), CAZymes (run_dbCAN), proteases (MEROPS),
// virulence (PHI-base / DFVF).
process EXTRAS {
  tag { meta.id }
  label 'extras'
  publishDir { "${params.outdir}/${meta.id}/13_extras" }, mode: 'copy'
  input:  tuple val(meta), path(proteins), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.extras.json"), emit: json
  script:
  """
  python3 ${projectDir}/bin/extras.py --sample "${meta.id}" --proteins ${proteins} \\
      --ont ${ont} --data-dir "${params.data_dir ?: ''}" --out ${meta.id}.extras.json \\
      || echo '{"sample":"${meta.id}","stage":"extras"}' > ${meta.id}.extras.json
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"extras\"}' > ${meta.id}.extras.json"
}
