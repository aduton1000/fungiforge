// Stage 03 — polishing. Medaka (ONT) always; when Illumina present and
// --hybrid != off, add short-read Polypolish + POLCA (hybrid branch). A
// per-base polishing-confidence flag rides downstream to the resistance caller.
process POLISH {
  tag { meta.id }
  label 'polish'
  publishDir { "${params.outdir}/${meta.id}/03_polish" }, mode: 'copy'
  input:  tuple val(meta), path(assembly), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.polished.fasta"), emit: assembly
          tuple val(meta), path("${meta.id}.polish.json"),    emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def do_hybrid = params.hybrid == 'on' || (params.hybrid == 'auto' && has_illumina)
  """
  mkdir -p medaka
  medaka_consensus -i ${ont} -d ${assembly} -o medaka -t ${task.cpus} -m ${params.dorado_model} || true
  [ -s medaka/consensus.fasta ] || cp ${assembly} medaka/consensus.fasta
  if [ "${do_hybrid}" = "true" ]; then
    bwa-mem2 index medaka/consensus.fasta 2>/dev/null || bwa index medaka/consensus.fasta
    ( bwa-mem2 mem -t ${task.cpus} medaka/consensus.fasta ${r1} > r1.sam && \\
      bwa-mem2 mem -t ${task.cpus} medaka/consensus.fasta ${r2} > r2.sam && \\
      polypolish polish medaka/consensus.fasta r1.sam r2.sam > ${meta.id}.polished.fasta ) || cp medaka/consensus.fasta ${meta.id}.polished.fasta
    POLISH_MODE=hybrid
  else
    cp medaka/consensus.fasta ${meta.id}.polished.fasta
    POLISH_MODE=ont_only
  fi
  python3 -c "import json;json.dump({'sample':'${meta.id}','stage':'polish','mode':'\${POLISH_MODE}','resistance_confidence':('high' if '\${POLISH_MODE}'=='hybrid' else 'provisional_ont_only')},open('${meta.id}.polish.json','w'),indent=2)"
  """
  stub:
  "touch ${meta.id}.polished.fasta ${meta.id}.polish.json"
}
