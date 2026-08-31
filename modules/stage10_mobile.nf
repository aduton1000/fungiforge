// Stage 10 — mobile & repeat elements (fungal sense). TE landscape from the
// Stage 06 library (LTR Gypsy/Copia, TIR, LINE, Helitron) + RIP index; mycovirus /
// endogenous-viral-element scan (geNomad + RVDB RdRp, DNA-only-limited); mito
// mobile introns / homing endonucleases; optional HGT alien-index.
process MOBILE {
  tag { meta.id }
  label 'mge'
  publishDir { "${params.outdir}/${meta.id}/10_mobile" }, mode: 'copy'
  input:  tuple val(meta), path(telib), path(mito)
  output: tuple val(meta), path("${meta.id}.mobile.json"), emit: json
  script:
  """
  python3 ${projectDir}/bin/te_summary.py --telib ${telib} --out te_summary.json || echo '{}' > te_summary.json
  if [ "${params.skip_mge}" != "true" ] && [ -n "${params.data_dir ?: ''}" ]; then
    genomad end-to-end --cleanup ${mito} genomad_out "${params.data_dir}/genomad_db" 2>/dev/null || true
  fi
  python3 ${projectDir}/bin/mobile_merge.py --sample "${meta.id}" --te te_summary.json \\
      --genomad genomad_out --out ${meta.id}.mobile.json || echo '{"sample":"${meta.id}","stage":"mobile"}' > ${meta.id}.mobile.json
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"mobile\"}' > ${meta.id}.mobile.json"
}
