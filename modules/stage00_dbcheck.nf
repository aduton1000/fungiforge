// Stage 00 (run-level) — verify the reference databases BEFORE any compute is spent.
// One task per run. Checks every database the enabled stages need (completion markers
// written by bin/fetch_references.sh + key files), writes db_manifest.json (embedded in
// provenance.json), and fails the run early with a single clear list of what is missing.
// --allow_missing_db turns the failure into a warning (dependent stages then record `skipped`).
process DB_CHECK {
  tag 'databases'
  label 'report'
  publishDir { "${params.outdir}/pipeline_info" }, mode: 'copy'
  output: path('db_manifest.json'), emit: manifest
  script:
  // which databases are required depends on which stages are enabled
  def req = ['unite', 'busco', 'funannotate', 'fungamr'] +
            (params.skip_decontam ? [] : ['kraken2']) +
            (params.skip_bgc      ? [] : ['antismash']) +
            (params.genome_id     ? ['refseq_fungi'] : [])
  def overrides = (params.funannotate_db ? "--path funannotate=${params.funannotate_db} " : '') +
                  (params.antismash_db   ? "--path antismash=${params.antismash_db} "     : '')
  if (params.data_dir)
    """
    check_databases.py --data-dir "${params.data_dir}" --out db_manifest.json \\
        --require ${req.join(',')} ${overrides} ${params.allow_missing_db ? '--allow-missing' : ''}
    """
  else
    """
    printf '{"data_dir": null, "ok": false, "manifest_tsv": false, "databases": {}, "optional": {}, "problems": ["--data_dir not set"]}\\n' > db_manifest.json
    echo "[db_check] --data_dir not set: identification, decontamination, annotation, resistance and BGC stages need the reference databases (bin/fetch_references.sh)" >&2
    ${params.allow_missing_db ? 'echo "[db_check] continuing because --allow_missing_db is set (dependent stages will record skipped)" >&2' : 'exit 1'}
    """
  stub:
  """
  printf '{"data_dir": "stub", "ok": true, "manifest_tsv": false, "databases": {}, "optional": {}, "problems": [], "stub": true}\\n' > db_manifest.json
  """
}
