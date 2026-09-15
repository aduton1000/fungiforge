// Stage 08 — identification & taxonomy. ITSx (ITS1/5.8S/ITS2) + barrnap LSU/SSU +
// secondary markers -> UNITE/RefSeq classification; genome-level sourmash gather +
// skani vs type strains; MLST (Candida/Aspergillus, incl. C. auris clade). GCPSR
// multi-locus concordance -> a per-isolate species call with confidence.
process IDENTIFY {
  tag { meta.id }
  label 'identify'
  publishDir { "${params.outdir}/${meta.id}/08_identify" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear)
  output: tuple val(meta), path("${meta.id}.species.txt"),  emit: species
          tuple val(meta), path("${meta.id}.markers.fasta"),emit: markers
          tuple val(meta), path("${meta.id}.identify.json"),emit: json
  script:
  """
  source "${projectDir}/bin/ff_status.sh"; ff_init identify "${meta.id}" ${meta.id}.identify.json
  ff_version barrnap -- barrnap --version
  ff_version ITSx -- ITSx --version
  ff_version vsearch -- vsearch --version
  ff_version sourmash -- sourmash --version
  ff_version blastn -- blastn -version
  ff_version mlst -- mlst --version
  # rRNA operon -> short region so ITSx doesn't hit the HMMER >100 kb limit on chromosomes
  ff_run barrnap -- bash -c "barrnap --kingdom fun --threads ${task.cpus} ${nuclear} > rrna.gff 2>barrnap.log"
  ff_run extract_rrna -- python3 ${projectDir}/bin/extract_rrna_region.py --assembly ${nuclear} --gff rrna.gff --out rrna_region.fasta
  ff_run ITSx -- bash -c "ITSx -i rrna_region.fasta -o itsx --cpu ${task.cpus} --preserve T --save_regions all 2>itsx.log"
  for f in itsx.full.fasta itsx.ITS1.fasta itsx.ITS2.fasta; do [ -f "\$f" ] && cat "\$f"; done > its.fasta
  cat its.fasta rrna_region.fasta > ${meta.id}.markers.fasta
  # PRIMARY id: ITS -> UNITE
  UNITE=\$(ls "${params.data_dir}"/unite/sh_general_release*dynamic*.fasta 2>/dev/null | grep -v _dev | head -1)
  if [ -s its.fasta ] && [ -n "\$UNITE" ]; then
    ff_run vsearch -- bash -c "vsearch --usearch_global its.fasta --db '\$UNITE' --id 0.90 --maxaccepts 10 --top_hits_only --blast6out unite.b6 --threads ${task.cpus} 2>vsearch.log"
  elif [ ! -s its.fasta ]; then ff_skip vsearch "no ITS region found in the assembly"
  else ff_skip vsearch "UNITE database not found under --data_dir"; fi
  # SECONDARY (optional) genome-level sourmash gather
  if [ "${params.genome_id}" = "true" ] && [ -n "${params.data_dir ?: ''}" ]; then
    ff_run sourmash_sketch --optional -- bash -c "sourmash sketch dna -p k=31,scaled=1000 ${nuclear} -o ${meta.id}.sig 2>sourmash.log"
    if [ "\$FF_RC" -eq 0 ]; then
      ff_run sourmash_gather --optional -- bash -c "sourmash gather ${meta.id}.sig '${params.data_dir}/refseq_fungi'/*.zip -o gather.csv 2>>sourmash.log"
    fi
  else
    ff_skip sourmash "genome-level identification disabled (--genome_id false) or no --data_dir"
  fi
  ff_run id_classify -- python3 ${projectDir}/bin/id_classify.py --sample "${meta.id}" --its its.fasta \\
      --unite-b6 unite.b6 --gather gather.csv \\
      --out-species ${meta.id}.species.txt --out-json ${meta.id}.identify.json
  ff_finalize
  """
  stub:
  "printf 'unknown\\t0\\n' > ${meta.id}.species.txt; touch ${meta.id}.markers.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"identify\"}' > ${meta.id}.identify.json"
}
