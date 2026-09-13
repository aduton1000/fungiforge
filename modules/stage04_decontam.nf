// Stage 04 — decontamination + organelle split. Kraken2 classifies every contig;
// bacterial/archaeal/viral/human contigs are dropped and the domain composition +
// a fungal/non_fungal verdict are recorded (surfaced in master_fungi.tsv); the mitochondrial genome is
// separated (oatk/GetOrganelle) because its mobile introns matter in Stage 10.
process DECONTAM {
  tag { meta.id }
  label 'decontam'
  publishDir { "${params.outdir}/${meta.id}/04_decontam" }, mode: 'copy'
  input:  tuple val(meta), path(assembly)
  output: tuple val(meta), path("${meta.id}.nuclear.fasta"), emit: nuclear
          tuple val(meta), path("${meta.id}.mito.fasta"),    emit: mito
          tuple val(meta), path("${meta.id}.decontam.json"), emit: json
  script:
  """
  if [ "${params.skip_decontam}" = "true" ] || [ -z "${params.data_dir ?: ''}" ] || [ ! -f "${params.data_dir ?: ''}/kraken2/hash.k2d" ]; then
    cp ${assembly} ${meta.id}.nuclear.fasta
    python3 -c "import json;json.dump({'sample':'${meta.id}','stage':'decontam','nuclear':'${meta.id}.nuclear.fasta','mito':'${meta.id}.mito.fasta','verdict':'not_run','note':'decontamination skipped (--skip_decontam, or no kraken2 DB under --data_dir)'},open('${meta.id}.decontam.json','w'),indent=2)"
  else
    # Kraken2 per contig -> drop Bacteria/Archaea/Viruses/human; keep fungal + unclassified.
    # A non-fungal isolate (verdict non_fungal/human) is kept WHOLE and flagged — there is
    # nothing to decontaminate, the sample is simply not a fungus (see kraken_taxonomy.py).
    kraken2 --db "${params.data_dir}/kraken2" --threads ${task.cpus} ${assembly} \\
        --output k2.out --report k2.report 2>k2.log
    python3 ${projectDir}/bin/kraken_taxonomy.py classify-contigs k2.report k2.out ${assembly} \\
        --sample "${meta.id}" --mito ${meta.id}.mito.fasta \\
        --out-fasta ${meta.id}.nuclear.fasta --out-json ${meta.id}.decontam.json
  fi
  # organelle: extract mito (placeholder — oatk/GetOrganelle wired in milestone 3)
  touch ${meta.id}.mito.fasta
  """
  stub:
  "touch ${meta.id}.nuclear.fasta ${meta.id}.mito.fasta ${meta.id}.decontam.json"
}
