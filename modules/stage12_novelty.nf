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
  source "${projectDir}/bin/ff_status.sh"; ff_init novelty "${meta.id}" ${meta.id}.novelty.json
  ff_version skani -- skani --version
  # genome-ANI novelty needs reference GENOMES (FASTA) — use skani if a genome set is staged;
  # otherwise the ITS-distance signal from Stage 08 is used alone and the skip is recorded.
  if [ "${params.skip_novelty}" = "true" ]; then ff_skip skani "disabled (--skip_novelty)"
  elif [ ! -d "${params.data_dir ?: ''}/refseq_fungi_genomes" ]; then ff_skip skani "reference genome set not staged under --data_dir (refseq_fungi_genomes)"
  else
    ff_run skani --optional -- bash -c "skani dist -q ${nuclear} -r ${params.data_dir}/refseq_fungi_genomes/*.f* -o skani.tsv 2>skani.log"
  fi
  ff_run novelty_call -- python3 ${projectDir}/bin/novelty_call.py --sample "${meta.id}" --skani skani.tsv \\
      --identify-json ${identify_json} --out ${meta.id}.novelty.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"novelty\"}' > ${meta.id}.novelty.json"
}
