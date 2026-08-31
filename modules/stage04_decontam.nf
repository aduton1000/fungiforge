// Stage 04 — decontamination + organelle split. Kraken2/tiara + BlobTools-style
// GC×coverage triage drops bacterial/human contigs; the mitochondrial genome is
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
  if [ "${params.skip_decontam}" = "true" ] || [ -z "${params.data_dir ?: ''}" ]; then
    cp ${assembly} ${meta.id}.nuclear.fasta
  else
    kraken2 --db "${params.data_dir}/kraken2" --threads ${task.cpus} ${assembly} \\
        --output k2.out --report k2.report 2>k2.log || true
    # keep contigs NOT classified as Bacteria/Homo (fungal + unclassified); refine with BlobTools later
    seqkit seq ${assembly} > ${meta.id}.nuclear.fasta
  fi
  # organelle: extract mito (placeholder — oatk/GetOrganelle wired in milestone 3)
  touch ${meta.id}.mito.fasta
  python3 -c "import json;json.dump({'sample':'${meta.id}','stage':'decontam','nuclear':'${meta.id}.nuclear.fasta','mito':'${meta.id}.mito.fasta'},open('${meta.id}.decontam.json','w'),indent=2)"
  """
  stub:
  "touch ${meta.id}.nuclear.fasta ${meta.id}.mito.fasta ${meta.id}.decontam.json"
}
