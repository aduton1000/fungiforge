// Stage 10 — mobile elements (W2.8). TE landscape: RepeatMasker's summary table (te_percent and
// the class breakdown) plus the RepeatModeler library (families by superfamily). geNomad on the
// whole assembly (nuclear + mitochondrial contigs): viruses, proviruses, plasmids (needs
// <data_dir>/genomad_db, fetch_references.sh genomad). Mycovirus / endogenous-viral-element
// screen: DIAMOND blastx of every contig against RVDB-prot (<data_dir>/rvdb, database built in
// the task), hits classified by virus family name; retrotransposon-like hits are reported apart.
// Mitochondrial mobile introns: homing-endonuclease ORFs (LAGLIDADG / GIY-YIG / HNH Pfam profiles
// from the funannotate database) on the mitogenome. Each block is optional (skipped with a reason).
process MOBILE {
  tag { meta.id }
  label 'mge'
  publishDir { "${params.outdir}/${meta.id}/10_mobile" }, mode: 'copy', pattern: "*.{json,tsv}"
  input:  tuple val(meta), path(telib), path(tbl), path(nuclear), path(mito)
  output: tuple val(meta), path("${meta.id}.mobile.json"), emit: json
  script:
  def data = params.data_dir ?: ''
  """
  source ff_status.sh; ff_init mobile "${meta.id}" ${meta.id}.mobile.json --best-effort
  ff_version genomad -- genomad --version
  ff_version diamond -- diamond version
  ff_version hmmer -- bash -c "hmmsearch -h | grep -m1 HMMER"
  ff_run te_summary -- te_summary.py --telib ${telib} --out te_summary.json
  if [ -s ${tbl} ]; then ff_run repeatmasker_tbl -- repeatmasker_tbl.py --tbl ${tbl} --out te_landscape.json
  else ff_skip repeatmasker_tbl "no RepeatMasker table from stage 06"; fi
  cat ${nuclear} ${mito} > all_contigs.fa
  # geNomad on the whole assembly
  if [ "${params.genomad}" != "true" ]; then ff_skip genomad "disabled (--genomad false)"
  elif [ -z "${data}" ] || [ ! -d "${data}/genomad_db" ]; then ff_skip genomad "geNomad database not staged under --data_dir/genomad_db (bin/fetch_references.sh genomad)"
  else
    ff_run genomad --optional -- bash -c "genomad end-to-end --cleanup --threads ${task.cpus} --splits ${params.genomad_splits} all_contigs.fa genomad_out '${data}/genomad_db' > genomad.log 2>&1"
  fi
  # RVDB-prot mycovirus / EVE screen
  MYCO=""
  if [ -s "${data}/rvdb/rvdb.fasta.xz" ] || [ -s "${data}/rvdb/rvdb.fasta" ]; then
    ff_run rvdb_makedb --optional -- bash -o pipefail -c "if [ -s '${data}/rvdb/rvdb.fasta.xz' ]; then xz -dc '${data}/rvdb/rvdb.fasta.xz'; else cat '${data}/rvdb/rvdb.fasta'; fi | diamond makedb --db rvdb --quiet --in -"
    if [ "\$FF_RC" -eq 0 ]; then
      ff_run rvdb_blastx --optional -- bash -c "diamond blastx -q all_contigs.fa -d rvdb -o rvdb_hits.tsv -p ${task.cpus} -e ${params.rvdb_evalue} -k 5 --quiet ${params.rvdb_sensitivity} --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle"
      if [ "\$FF_RC" -eq 0 ]; then
        MITO_IDS=\$(grep '>' ${mito} | sed 's/>//; s/ .*//' | paste -sd, -)
        ff_run mycovirus_screen -- mycovirus_screen.py --sample "${meta.id}" --hits rvdb_hits.tsv --mito-contigs "\${MITO_IDS:-}" --out mycovirus.json
        MYCO="--mycovirus mycovirus.json"; cp rvdb_hits.tsv ${meta.id}.rvdb_hits.tsv
      fi
    fi
  else
    ff_skip rvdb "RVDB-prot not staged under --data_dir/rvdb (bin/fetch_references.sh rvdb)"
  fi
  # mitochondrial homing-endonuclease ORFs
  PFAM="${params.funannotate_db ?: data + '/funannotate'}/Pfam-A.hmm"
  if [ -s ${mito} ]; then
    ff_run mito_heg --optional -- mito_heg.py --sample "${meta.id}" --mito ${mito} --pfam-hmm "\$PFAM" --threads ${task.cpus} --workdir heg_work --out heg.json
    HEG=\$([ -s heg.json ] && echo "--heg heg.json" || echo "")
  else
    ff_skip mito_heg "no mitochondrial contigs"; HEG=""
  fi
  ff_run mobile_merge -- mobile_merge.py --sample "${meta.id}" --te te_summary.json --tbl te_landscape.json \\
      --genomad genomad_out \$MYCO \$HEG --out ${meta.id}.mobile.json
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"mobile\"}' > ${meta.id}.mobile.json"
}
