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
  ITSx -i ${nuclear} -o itsx --cpu ${task.cpus} --preserve T || true
  barrnap --kingdom euk --threads ${task.cpus} ${nuclear} > rrna.gff 2>/dev/null || true
  cat itsx*.fasta rrna*.fasta 2>/dev/null > ${meta.id}.markers.fasta || touch ${meta.id}.markers.fasta
  # genome-level: sourmash gather + skani vs reference (needs data_dir)
  if [ -n "${params.data_dir ?: ''}" ]; then
    sourmash sketch dna -p k=31,scaled=1000 ${nuclear} -o ${meta.id}.sig 2>/dev/null || true
    sourmash gather ${meta.id}.sig "${params.data_dir}/refseq_fungi"/*.zip -o gather.csv 2>/dev/null || true
  fi
  python3 ${projectDir}/bin/id_classify.py --sample "${meta.id}" --markers ${meta.id}.markers.fasta \\
      --gather gather.csv --data-dir "${params.data_dir ?: ''}" \\
      --out-species ${meta.id}.species.txt --out-json ${meta.id}.identify.json || \\
      { echo "unknown\t0" > ${meta.id}.species.txt; echo '{"sample":"${meta.id}","stage":"identify","species":"unknown"}' > ${meta.id}.identify.json; }
  """
  stub:
  "printf 'unknown\\t0\\n' > ${meta.id}.species.txt; touch ${meta.id}.markers.fasta ${meta.id}.identify.json"
}
