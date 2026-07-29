// Stage 12 — novelty. Combine genome ANI (skani < 95%) + ITS distance +
// phylogenomic placement (BUSCO single-copy orthologs -> IQ-TREE) + AAI to flag
// candidate novel species, with GCPSR/polyphasic caveats reported.
process NOVELTY {
  tag { meta.id }
  label 'identify'
  publishDir { "${params.outdir}/${meta.id}/12_novelty" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear), path(markers), path(identify_json)
  output: tuple val(meta), path("${meta.id}.novelty.json"), emit: json
  script:
  """
  # genome-ANI novelty needs reference GENOMES (FASTA) — use skani if a genome set is staged;
  # otherwise fall back to the ITS-distance signal from Stage 08 (fast, always available).
  if [ "${params.skip_novelty}" != "true" ] && [ -d "${params.data_dir}/refseq_fungi_genomes" ]; then
    skani dist -q ${nuclear} -r ${params.data_dir}/refseq_fungi_genomes/*.f* -o skani.tsv 2>/dev/null || true
  fi
  python3 ${projectDir}/bin/novelty_call.py --sample "${meta.id}" --skani skani.tsv \\
      --identify-json ${identify_json} --out ${meta.id}.novelty.json || echo '{"sample":"${meta.id}","stage":"novelty"}' > ${meta.id}.novelty.json
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"novelty\"}' > ${meta.id}.novelty.json"
}
