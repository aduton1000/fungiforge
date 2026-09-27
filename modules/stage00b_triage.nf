// Stage 00b — read-level triage (W2.2). Kraken2 on a subsample of the QC'd reads gives an
// early verdict (fungal | likely_fungal | mixed | non_fungal | human | not_run) with the same
// rules the contig-level verdict uses after assembly (bin/kraken_taxonomy.py). An isolate that
// is non_fungal or human at read level is stopped by the gate BEFORE assembly, so a bacterium
// on a fungal plate costs a few minutes of Kraken2, not hours of SPAdes/Flye and annotation.
// Illumina reads are used when present (hybrid/short-read), ONT otherwise. Skipped (verdict
// not_run) when no Kraken2 database is staged.
process READ_TRIAGE {
  tag { meta.id }
  label 'decontam'
  publishDir { "${params.outdir}/${meta.id}/01_readqc" }, mode: 'copy', pattern: "*.json"
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.triage.json"), emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def platform     = has_illumina ? 'illumina' : 'ont'
  def n            = params.triage_reads
  def sub_cmd      = has_illumina ? "seqkit head -n ${n} ${r1} > sub_1.fq && seqkit head -n ${n} ${r2} > sub_2.fq"
                                  : "seqkit head -n ${n} ${ont} > sub.fq"
  def total_src    = has_illumina ? r1 : ont
  def k2_in        = has_illumina ? "--paired sub_1.fq sub_2.fq" : "sub.fq"
  def sub_file     = has_illumina ? "sub_1.fq" : "sub.fq"
  """
  source ff_status.sh; ff_init triage "${meta.id}" ${meta.id}.triage.json
  ff_version kraken2 -- kraken2 --version
  ff_version seqkit -- seqkit version
  if [ -z "${params.data_dir ?: ''}" ] || [ ! -f "${params.data_dir ?: ''}/kraken2/hash.k2d" ]; then
    ff_skip kraken2 "no kraken2 database under --data_dir: read triage not run"
    ff_run read_triage -- read_triage.py --sample "${meta.id}" --report /dev/null --n-reads 0 --platform ${platform} --out ${meta.id}.triage.json
  else
    ff_run subsample -- bash -o pipefail -c "${sub_cmd}"
    TOTAL=\$(seqkit stats -T ${total_src} | awk 'NR==2{print \$4}')
    ff_run kraken2 -- bash -c "kraken2 --db '${params.data_dir}/kraken2' --threads ${task.cpus} ${k2_in} --output /dev/null --report k2.report 2>k2.log"
    NSUB=\$(grep -c '^+\$' ${sub_file} || true)
    ff_run read_triage -- read_triage.py --sample "${meta.id}" --report k2.report --n-reads "\${NSUB:-0}" \\
        --total-reads "\${TOTAL:-0}" --platform ${platform} --out ${meta.id}.triage.json
  fi
  ff_finalize
  """
  stub:
  // --stub_read_verdicts 'ID=non_fungal,ID2=human' lets the DAG tests exercise the read-level gate
  def forced = (params.stub_read_verdicts ?: '').toString().split(',').collect { kv -> kv.split('=') }.findAll { kv -> kv.size() == 2 }.collectEntries { kv -> [(kv[0].trim()): kv[1].trim()] }
  def verdict = forced[meta.id] ?: 'likely_fungal'
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"triage\",\"verdict\":\"${verdict}\"}' > ${meta.id}.triage.json"
}
