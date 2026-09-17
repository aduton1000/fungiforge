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
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"   # node /tmp is a shared tmpfs; keep scratch on disk
  source ff_status.sh; ff_init assembly_qc "${meta.id}" ${meta.id}.assemblyqc.json
  ff_version busco -- busco --version
    ff_version compleasm -- compleasm --version
  # Contiguity statistics come from assembly_qc.py (no QUAST in the images; the reference
  # benchmark lives in bin/benchmark_assembly.py, run by the validation suite, not per isolate).
  COMPLEASM_OK=0
  if command -v compleasm >/dev/null 2>&1 && [ -n "${params.data_dir ?: ''}" ]; then
    ff_run compleasm --optional -- compleasm run -a ${nuclear} -o compleasm -l ${lineage} -L "${params.data_dir}/busco" -t ${task.cpus}
    if [ "\$FF_RC" -eq 0 ]; then COMPLEASM_OK=1; fi
  else
    ff_skip compleasm "compleasm not available in this image; using BUSCO"
  fi
  if [ "\$COMPLEASM_OK" -eq 0 ]; then
    ff_run busco -- busco -i ${nuclear} -o busco -l ${lineage} -m genome -c ${task.cpus} \\
        ${ params.data_dir ? "--download_path ${params.data_dir}/busco --offline" : "" }
  fi
  ff_run assembly_qc -- assembly_qc.py --sample "${meta.id}" --nuclear ${nuclear} \\
      --lineage ${lineage} --compleasm-dir compleasm --busco-dir busco \\
      --out ${meta.id}.assemblyqc.json
  ff_finalize
  """
  stub:
  // --stub_qc_fail 'ID,ID2' lets the DAG tests exercise the gate
  def failing = (params.stub_qc_fail ?: '').toString().split(',').collect { id -> id.trim() }
  def qc = failing.contains(meta.id) ? 'false' : 'true'
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"assembly_qc\",\"qc_pass\":${qc}}' > ${meta.id}.assemblyqc.json"
}
