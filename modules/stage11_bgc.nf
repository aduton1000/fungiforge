// Stage 11 — biosynthetic gene clusters (fungiSMASH / antiSMASH 8 fungal mode).
// NRPS/PKS/terpene/RiPP/hybrid clusters; downstream BiG-SCAPE networks + MIBiG
// (known mycotoxin/antibiotic vs novel) run at the comparative layer across isolates.
process BGC {
  tag { meta.id }
  label 'bgc'
  publishDir { "${params.outdir}/${meta.id}/11_bgc" }, mode: 'copy'
  input:  tuple val(meta), path(gbk)
  output: tuple val(meta), path("${meta.id}.bgc.json"),      emit: json
          tuple val(meta), path("${meta.id}.regions.gbk"),   emit: regions, optional: true
  script:
  """
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"   # node /tmp is a shared tmpfs; keep scratch on disk
  # BGC is best-effort: antiSMASH needs a GBK with annotated CDS (genefinding none). An
  # empty/geneless GBK (e.g. failed annotation) makes antiSMASH exit non-zero — that must
  # NOT kill a multi-hour run, so we guard the input and force a clean exit 0.
  set +e
  mkdir -p as
  RUN_BGC=0
  if [ "${params.skip_bgc}" != "true" ] && [ -s ${gbk} ] && grep -q "^LOCUS" ${gbk} 2>/dev/null; then
    if grep -qE "^     CDS |/translation=" ${gbk} 2>/dev/null; then RUN_BGC=1; fi
  fi
  STATUS=skipped
  if [ "\$RUN_BGC" = "1" ]; then
    antismash --taxon fungi --output-dir as --genefinding-tool none --cpus ${task.cpus} \\
        --databases "${params.antismash_db ?: params.data_dir + '/antismash'}" ${gbk} > as.log 2>&1
    RC=\$?
    if [ \$RC -eq 0 ]; then STATUS=ok; else STATUS=failed; echo "antismash exited \$RC — see as.log" >&2; tail -20 as.log >&2; fi
    cat as/*.region*.gbk > ${meta.id}.regions.gbk 2>/dev/null || true
  else
    echo "BGC skipped: GBK empty or has no annotated CDS (skip_bgc=${params.skip_bgc})." >&2
  fi
  # Best-effort stage, but never a silent zero: a failed/skipped antiSMASH run is recorded as
  # status failed/skipped with n_clusters null (-> NA in master_fungi.tsv), not as "0 BGCs".
  python3 ${projectDir}/bin/bgc_summary.py --sample "${meta.id}" --as-dir as --status "\$STATUS" \\
      --log as.log --out ${meta.id}.bgc.json \\
    || printf '{"sample":"%s","stage":"bgc","status":"failed","n_clusters":null,"clusters":[],"note":"bgc_summary failed"}\\n' "${meta.id}" > ${meta.id}.bgc.json
  exit 0
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"bgc\"}' > ${meta.id}.bgc.json"
}
