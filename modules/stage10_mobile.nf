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
  source "${projectDir}/bin/ff_status.sh"; ff_init mobile "${meta.id}" ${meta.id}.mobile.json
  ff_version genomad -- genomad --version
  ff_run te_summary -- python3 ${projectDir}/bin/te_summary.py --telib ${telib} --out te_summary.json
  if [ "${params.skip_mge}" = "true" ]; then ff_skip genomad "disabled (--skip_mge)"
  elif [ -z "${params.data_dir ?: ''}" ] || [ ! -d "${params.data_dir ?: ''}/genomad_db" ]; then ff_skip genomad "geNomad database not staged under --data_dir (genomad_db)"
  else
    ff_run genomad --optional -- bash -c "genomad end-to-end --cleanup ${mito} genomad_out '${params.data_dir}/genomad_db' > genomad.log 2>&1"
  fi
  ff_run mobile_merge -- python3 ${projectDir}/bin/mobile_merge.py --sample "${meta.id}" --te te_summary.json \\
      --genomad genomad_out --out ${meta.id}.mobile.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"mobile\"}' > ${meta.id}.mobile.json"
}
