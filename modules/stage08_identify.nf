// Stage 08 — identification & taxonomy (W2.3). ITSx (ITS1/5.8S/ITS2) against UNITE is the
// primary barcode; the secondary loci (CaM, BenA, TEF1, RPB2 located by tblastn of bundled
// reference proteins, LSU D1/D2 from barrnap's 28S) are searched against per-locus type-material
// reference sets (<data_dir>/markers, bin/fetch_references.sh markers); optional genome-level
// sourmash gather; MLST against the fungal PubMLST schemes (<data_dir>/mlst) when staged.
// bin/id_classify.py applies the concordance rule (species needs ITS + one agreeing secondary
// line; discordance -> genus + flag) and chooses the species-aware BUSCO lineage for stage 08b.
// Absent databases skip their step (status partial/skipped), never fail the isolate.
process IDENTIFY {
  tag { meta.id }
  label 'identify'
  publishDir { "${params.outdir}/${meta.id}/08_identify" }, mode: 'copy'
  input:  tuple val(meta), path(nuclear)
  output: tuple val(meta), path("${meta.id}.species.txt"),       emit: species
          tuple val(meta), path("${meta.id}.markers.fasta"),     emit: markers
          tuple val(meta), path("${meta.id}.identify.json"),     emit: json
          tuple val(meta), path("${meta.id}.busco_lineage.txt"), emit: lineage
  script:
  def data = params.data_dir ?: ''
  """
  source ff_status.sh; ff_init identify "${meta.id}" ${meta.id}.identify.json
  ff_version barrnap -- barrnap --version
  ff_version ITSx -- ITSx --version
  ff_version vsearch -- vsearch --version
  ff_version sourmash -- sourmash --version
  ff_version blastn -- blastn -version
  ff_version mlst -- mlst --version
  # rRNA operon -> short region so ITSx doesn't hit the HMMER >100 kb limit on chromosomes
  ff_run barrnap -- bash -c "barrnap --kingdom fun --threads ${task.cpus} ${nuclear} > rrna.gff 2>barrnap.log"
  ff_run extract_rrna -- extract_rrna_region.py --assembly ${nuclear} --gff rrna.gff --out rrna_region.fasta
  ff_run ITSx -- bash -c "ITSx -i rrna_region.fasta -o itsx --cpu ${task.cpus} --preserve T --save_regions all 2>itsx.log"
  for f in itsx.full.fasta itsx.ITS1.fasta itsx.ITS2.fasta; do [ -f "\$f" ] && cat "\$f"; done > its.fasta
  # PRIMARY id: ITS -> UNITE
  UNITE=\$(ls "${data}"/unite/sh_general_release*dynamic*.fasta 2>/dev/null | grep -v _dev | head -1)
  if [ -s its.fasta ] && [ -n "\$UNITE" ]; then
    ff_run vsearch -- bash -c "vsearch --usearch_global its.fasta --db '\$UNITE' --id 0.90 --maxaccepts 10 --top_hits_only --blast6out unite.b6 --threads ${task.cpus} 2>vsearch.log"
  elif [ ! -s its.fasta ]; then ff_skip vsearch "no ITS region found in the assembly"
  else ff_skip vsearch "UNITE database not found under --data_dir"; fi
  # SECONDARY loci: locate with tblastn (bundled reference proteins), extract, search the type-material sets
  ff_run makeblastdb -- bash -c "makeblastdb -in ${nuclear} -dbtype nucl -out asm_db >makeblastdb.log 2>&1"
  ff_run tblastn -- bash -c "tblastn -query '${params.marker_proteins}' -db asm_db -evalue 1e-10 -max_target_seqs 20 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > tblastn.tsv 2>tblastn.log"
  ff_run extract_markers -- extract_markers.py --sample "${meta.id}" --assembly ${nuclear} --tblastn tblastn.tsv --rrna-gff rrna.gff \\
      --flank ${params.marker_flank} --out-prefix ${meta.id}.marker --json markers_extracted.json
  LOCUS_ARGS=""
  for locus in CaM BenA TEF1 RPB2 LSU; do
    q=${meta.id}.marker.\$locus.fasta; ref="${data}/markers/\$locus.fasta"
    if [ -s "\$q" ] && [ -s "\$ref" ]; then
      ff_run blastn_\$locus --optional -- bash -c "makeblastdb -in '\$ref' -dbtype nucl -out ref_\$locus >>makeblastdb.log 2>&1 && blastn -query \$q -db ref_\$locus -task megablast -evalue 1e-20 -max_target_seqs 50 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle' > \$locus.b6 2>>blastn.log"
      [ "\$FF_RC" -eq 0 ] && LOCUS_ARGS="\$LOCUS_ARGS --locus-b6 \$locus=\$locus.b6"
    elif [ ! -s "\$q" ]; then ff_skip blastn_\$locus "\$locus not found in the assembly"
    else ff_skip blastn_\$locus "no \$locus reference set under --data_dir/markers"; fi
  done
  cat its.fasta rrna_region.fasta ${meta.id}.marker.*.fasta > ${meta.id}.markers.fasta
  # MLST: fungal PubMLST schemes staged by fetch_references.sh mlst; the BLAST db is built here
  if ls "${data}"/mlst/pubmlst/*/*.tfa >/dev/null 2>&1; then
    ff_run mlst --optional -- mlst_fungal.sh "${data}/mlst/pubmlst" ${nuclear} ${task.cpus} mlst.tsv
    MLST_ARG="--mlst mlst.tsv"; [ "\$FF_RC" -eq 0 ] || MLST_ARG=""
  else
    ff_skip mlst "no PubMLST schemes under --data_dir/mlst (bin/fetch_references.sh mlst)"; MLST_ARG=""
  fi
  # OPTIONAL genome-level sourmash gather
  if [ "${params.genome_id}" = "true" ] && [ -n "${data}" ]; then
    ff_run sourmash_sketch --optional -- bash -c "sourmash sketch dna -p k=31,scaled=1000 ${nuclear} -o ${meta.id}.sig 2>sourmash.log"
    if [ "\$FF_RC" -eq 0 ]; then
      ff_run sourmash_gather --optional -- bash -c "sourmash gather ${meta.id}.sig '${data}/refseq_fungi'/*.zip -o gather.csv 2>>sourmash.log"
    fi
  else
    ff_skip sourmash "genome-level identification disabled (--genome_id false) or no --data_dir"
  fi
  ff_run id_classify -- id_classify.py --sample "${meta.id}" --its its.fasta \\
      --unite-b6 unite.b6 --gather gather.csv --markers-json markers_extracted.json \$LOCUS_ARGS \$MLST_ARG \\
      --lineage-map "${params.busco_lineage_map}" \\
      --out-species ${meta.id}.species.txt --out-json ${meta.id}.identify.json --out-lineage ${meta.id}.busco_lineage.txt
  ff_finalize
  """
  stub:
  """
  printf 'unknown\\t0\\n' > ${meta.id}.species.txt; touch ${meta.id}.markers.fasta
  echo '{"sample":"${meta.id}","stage":"identify","species":"unknown","busco_lineage":"${params.stub_busco_lineage ?: 'fungi_odb10'}"}' > ${meta.id}.identify.json
  echo '${params.stub_busco_lineage ?: 'fungi_odb10'}' > ${meta.id}.busco_lineage.txt
  """
}
