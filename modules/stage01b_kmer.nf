// Stage 01b — k-mer profile (W2.2): KMC k-mer spectrum + GenomeScope2 (ploidy 1 and 2) on
// the QC'd reads -> estimated genome size, heterozygosity, ploidy hint, k-mer coverage and
// read coverage (bin/kmer_profile.py). Trimmed Illumina reads are profiled when present
// (accurate k-mers), else the filtered ONT reads (noisier; models may not converge — then the
// values are null with a note, never invented). The tools are optional steps: a failed fit
// leaves the stage `partial` and the run continues; the profile JSON is always written.
process KMER_PROFILE {
  tag { meta.id }
  label 'readqc'
  publishDir { "${params.outdir}/${meta.id}/01_readqc" }, mode: 'copy', pattern: "*.json"
  input:  tuple val(meta), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.kmer.json"), emit: json
  script:
  def has_illumina = r1.name != 'NO_R1'
  def source       = has_illumina ? 'illumina_trimmed' : 'ont_filtered'
  def files        = has_illumina ? "${r1} ${r2}" : "${ont}"
  def k            = params.kmer_k
  def mem_gb       = Math.max(4, (task.memory.toGiga() - 2) as int)
  """
  source ff_status.sh; ff_init kmer "${meta.id}" ${meta.id}.kmer.json
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR" kmc_tmp
  ff_version kmc -- bash -c "kmc 2>&1 | head -1"
  ff_version genomescope2 -- genomescope2 --version
  ff_version seqkit -- seqkit version
  printf '%s\\n' ${files} > kmc_in.txt
  # read totals for the coverage estimate (sum of bases, length-weighted mean read length)
  seqkit stats -T ${files} | awk 'NR>1 {b+=\$5; n+=\$4} END {printf "%d\\t%.1f\\n", b, (n ? b/n : 0)}' > read_totals.tsv
  BASES=\$(cut -f1 read_totals.tsv); MEANLEN=\$(cut -f2 read_totals.tsv)
  ff_run kmc --optional -- bash -c "kmc -k${k} -t${task.cpus} -m${mem_gb} -ci1 -cs100000 -fq @kmc_in.txt kmc_db kmc_tmp > kmc.log 2>&1"
  if [ "\$FF_RC" -eq 0 ]; then
    ff_run kmc_hist --optional -- bash -c "kmc_tools transform kmc_db histogram kmc.hist -cx100000 > kmc_hist.log 2>&1"
    ff_run genomescope_p1 --optional -- bash -c "genomescope2 -i kmc.hist -o gs_p1 -k ${k} -p 1 > gs_p1.log 2>&1"
    ff_run genomescope_p2 --optional -- bash -c "genomescope2 -i kmc.hist -o gs_p2 -k ${k} -p 2 > gs_p2.log 2>&1"
  else
    ff_skip genomescope2 "k-mer counting failed; no spectrum to model"
  fi
  rm -rf kmc_tmp kmc_db.kmc_pre kmc_db.kmc_suf
  ff_run kmer_profile -- kmer_profile.py --sample "${meta.id}" --hist kmc.hist --gs-p1 gs_p1 --gs-p2 gs_p2 --k ${k} \\
      --read-bases "\${BASES:-0}" --mean-read-len "\${MEANLEN:-0}" --source ${source} --out ${meta.id}.kmer.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"kmer\",\"genome_size_est\":29000000,\"heterozygosity_pct\":0.1,\"ploidy_hint\":\"haploid\",\"coverage_from_kmers\":60.0}' > ${meta.id}.kmer.json"
}
