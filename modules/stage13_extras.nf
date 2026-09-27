// Stage 13 — eukaryote-specific extras (W2.6). Secretome (SignalP 6 in the site-built image) and
// effectors (EffectorP 3), CAZymes (dbCAN HMMs with hmmsearch), virulence (PHI-base, diamond),
// mating type (Pfam MATalpha_HMGbox / HMG_box in the MAT locus context of the GenBank), ploidy
// from the read BAM of stage 09 (nQuire). Each block is optional: a missing tool or database is
// recorded in `tools`/`skipped` and the stage is `partial`, never failed (bin/extras.py).
process EXTRAS {
  tag { meta.id }
  label 'extras'
  publishDir { "${params.outdir}/${meta.id}/13_extras" }, mode: 'copy', pattern: "*.json"
  input:  tuple val(meta), path(proteins), path(gbk), path(bam), path(bai)
  output: tuple val(meta), path("${meta.id}.extras.json"), emit: json
  script:
  def bam_arg = bam.name.endsWith('.bam') ? "--bam ${bam}" : ''
  """
  source ff_status.sh; ff_init extras "${meta.id}" ${meta.id}.extras.json --best-effort
  ff_version python -- python3 --version
  ff_version hmmer -- bash -c "hmmsearch -h | grep -m1 HMMER"
  ff_version diamond -- diamond version
  ff_version nquire -- bash -c "nQuire 2>&1 | head -1"
  ff_version signalp6 -- bash -c "signalp6 --version 2>&1 | head -1"
  ff_version java -- bash -c "java -version 2>&1 | head -1"
  ff_run extras -- extras.py --sample "${meta.id}" --proteins ${proteins} --gbk ${gbk} ${bam_arg} \\
      --data-dir "${params.data_dir ?: ''}" --threads ${task.cpus} --workdir extras_work \
      --ploidy-min-sites ${params.ploidy_min_sites} --out ${meta.id}.extras.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"extras\",\"mating_type\":\"undetermined\",\"ploidy\":\"NA\"}' > ${meta.id}.extras.json"
}
