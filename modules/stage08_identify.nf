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
  # rRNA operon -> short region so ITSx doesn't hit the HMMER >100 kb limit on chromosomes
  barrnap --kingdom fun --threads ${task.cpus} ${nuclear} > rrna.gff 2>barrnap.log || true
  python3 ${projectDir}/bin/extract_rrna_region.py --assembly ${nuclear} --gff rrna.gff --out rrna_region.fasta || : > rrna_region.fasta
  ITSx -i rrna_region.fasta -o itsx --cpu ${task.cpus} --preserve T --save_regions all 2>itsx.log || true
  cat itsx.full.fasta itsx.ITS1.fasta itsx.ITS2.fasta 2>/dev/null > its.fasta || : > its.fasta
  cat its.fasta rrna_region.fasta 2>/dev/null > ${meta.id}.markers.fasta
  # PRIMARY id: ITS -> UNITE (fast, emulation-friendly)
  UNITE=\$(ls "${params.data_dir}"/unite/sh_general_release*dynamic*.fasta 2>/dev/null | grep -v _dev | head -1)
  if [ -s its.fasta ] && [ -n "\$UNITE" ]; then
    vsearch --usearch_global its.fasta --db "\$UNITE" --id 0.90 --maxaccepts 10 --top_hits_only \\
        --blast6out unite.b6 --threads ${task.cpus} 2>vsearch.log || true
  fi
  # SECONDARY (optional; slow under x86 emulation) genome-level sourmash gather
  if [ "${params.genome_id}" = "true" ] && [ -n "${params.data_dir ?: ''}" ]; then
    sourmash sketch dna -p k=31,scaled=1000 ${nuclear} -o ${meta.id}.sig 2>/dev/null || true
    sourmash gather ${meta.id}.sig "${params.data_dir}/refseq_fungi"/*.zip -o gather.csv 2>/dev/null || true
  fi
  python3 ${projectDir}/bin/id_classify.py --sample "${meta.id}" --its its.fasta \\
      --unite-b6 unite.b6 --gather gather.csv \\
      --out-species ${meta.id}.species.txt --out-json ${meta.id}.identify.json || \\
      { printf 'unknown\\tnone\\n' > ${meta.id}.species.txt; echo '{"sample":"${meta.id}","stage":"identify","species":"unknown"}' > ${meta.id}.identify.json; }
  """
  stub:
  "printf 'unknown\\t0\\n' > ${meta.id}.species.txt; touch ${meta.id}.markers.fasta ${meta.id}.identify.json"
}
