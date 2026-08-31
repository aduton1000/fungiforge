// Stage 05 — assembly QC + completeness. QUAST contiguity + compleasm/BUSCO
// against the order-specific ODB lineage (auto after ID; fungi_odb10 fallback).
// Emits a QC-pass flag (MIMAG-style) carried into the comparative layer.
process ASSEMBLY_QC {
  tag { meta.id }
  label 'busco'
  publishDir { "${params.outdir}/${meta.id}/05_assembly_qc" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear)
  output: tuple val(meta), path("${meta.id}.assemblyqc.json"), emit: json
  script:
  def lineage = params.busco_lineage == 'auto' ? 'fungi_odb10' : params.busco_lineage
  """
  quast.py ${nuclear} -o quast --threads ${task.cpus} --silent || true
  if command -v compleasm >/dev/null 2>&1 && [ -n "${params.data_dir ?: ''}" ]; then
    compleasm run -a ${nuclear} -o compleasm -l ${lineage} -L "${params.data_dir}/busco" -t ${task.cpus} || true
  else
    busco -i ${nuclear} -o busco -l ${lineage} -m genome -c ${task.cpus} \\
        ${ params.data_dir ? "--download_path ${params.data_dir}/busco --offline" : "" } || true
  fi
  python3 ${projectDir}/bin/assembly_qc.py --sample "${meta.id}" --nuclear ${nuclear} \\
      --lineage ${lineage} --compleasm-dir compleasm --busco-dir busco \\
      --out ${meta.id}.assemblyqc.json
  """
  stub:
  "touch ${meta.id}.assemblyqc.json"
}
