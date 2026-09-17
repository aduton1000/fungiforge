// Stage 12 — novelty (W2.9). skani ANI against the staged reference genome set
// (<data_dir>/refseq_fungi_genomes, bin/fetch_references.sh genomes: the reference genome of
// every species in the genera of fungiforge/resources/novelty_genera.txt) combined with the
// ITS identity and the secondary-locus concordance of stage 08: known_species (ANI >= 95 %),
// candidate_novel_species (90-95 %), candidate_novel_or_unrepresented (< 90 % or no aligned
// reference); ITS-distance fallback when no genome set is staged.
process NOVELTY {
  tag { meta.id }
  label 'identify'
  publishDir { "${params.outdir}/${meta.id}/12_novelty" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear), path(markers), path(identify_json)
  output: tuple val(meta), path("${meta.id}.novelty.json"), emit: json
  script:
  """
  source ff_status.sh; ff_init novelty "${meta.id}" ${meta.id}.novelty.json
  ff_version skani -- skani --version
  # genome-ANI novelty needs reference GENOMES (FASTA) — use skani if a genome set is staged;
  # otherwise the ITS-distance signal from Stage 08 is used alone and the skip is recorded.
  if [ "${params.skip_novelty}" = "true" ]; then ff_skip skani "disabled (--skip_novelty)"
  elif ! ls "${params.data_dir ?: ''}"/refseq_fungi_genomes/*.fna >/dev/null 2>&1; then ff_skip skani "reference genome set not staged under --data_dir/refseq_fungi_genomes (bin/fetch_references.sh genomes)"
  else
    ls "${params.data_dir}"/refseq_fungi_genomes/*.fna > refs.txt
    ff_run skani --optional -- bash -c "skani dist -q ${nuclear} --rl refs.txt -t ${task.cpus} --min-af ${params.novelty_min_af} -o skani.tsv 2>skani.log"
  fi
  MANIFEST="${params.data_dir ?: ''}/refseq_fungi_genomes/manifest.tsv"
  ff_run novelty_call -- novelty_call.py --sample "${meta.id}" --skani skani.tsv --min-af ${params.novelty_min_af} \\
      --identify-json ${identify_json} \$([ -s "\$MANIFEST" ] && echo "--manifest \$MANIFEST") --out ${meta.id}.novelty.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"novelty\"}' > ${meta.id}.novelty.json"
}
