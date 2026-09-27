// Stage 04 — decontamination + organelle split. Kraken2 classifies every contig;
// bacterial/archaeal/viral/human contigs are dropped and the domain composition +
// a fungal/non_fungal verdict are recorded (surfaced in master_fungi.tsv); the mitochondrial genome is
// separated from the assembly (W2.7: contigs carrying >= 2 core mitochondrial genes by tblastn of the
// bundled reference proteins, <= 250 kb, AT-rich) because its mobile introns matter in Stage 10 and
// stage 04b annotates it.
process DECONTAM {
  tag { meta.id }
  label 'decontam'
  publishDir { "${params.outdir}/${meta.id}/04_decontam" }, mode: 'copy'
  input:  tuple val(meta), path(assembly)
  output: tuple val(meta), path("${meta.id}.nuclear.fasta"), emit: nuclear
          tuple val(meta), path("${meta.id}.mito.fasta"),    emit: mito
          tuple val(meta), path("${meta.id}.decontam.json"), emit: json
  script:
  """
  source ff_status.sh; ff_init decontam "${meta.id}" ${meta.id}.decontam.json
  ff_version kraken2 -- kraken2 --version
  if [ "${params.skip_decontam}" = "true" ] || [ -z "${params.data_dir ?: ''}" ] || [ ! -f "${params.data_dir ?: ''}/kraken2/hash.k2d" ]; then
    ff_skip kraken2 "decontamination skipped (--skip_decontam, or no kraken2 DB under --data_dir)"
    cp ${assembly} ${meta.id}.nuclear.fasta
    python3 -c "import json;json.dump({'sample':'${meta.id}','stage':'decontam','nuclear':'${meta.id}.nuclear.fasta','mito':'${meta.id}.mito.fasta','verdict':'not_run'},open('${meta.id}.decontam.json','w'),indent=2)"
  else
    # Kraken2 per contig -> drop Bacteria/Archaea/Viruses/human; keep fungal + unclassified.
    # A non-fungal isolate (verdict non_fungal/human) is kept WHOLE and flagged — there is
    # nothing to decontaminate, the sample is simply not a fungus (see kraken_taxonomy.py).
    ff_run kraken2 -- bash -c "kraken2 --db '${params.data_dir}/kraken2' --threads ${task.cpus} ${assembly} --output k2.out --report k2.report 2>k2.log"
    ff_run classify_contigs -- kraken_taxonomy.py classify-contigs k2.report k2.out ${assembly} \\
        --sample "${meta.id}" --mito ${meta.id}.mito.fasta \\
        --out-fasta ${meta.id}.nuclear.fasta --out-json ${meta.id}.decontam.json
  fi
  # organelle split (W2.7): the mitochondrial contigs leave the nuclear set
  ff_run makeblastdb -- bash -c "makeblastdb -in ${meta.id}.nuclear.fasta -dbtype nucl -out asm_db > makeblastdb.log 2>&1"
  ff_run tblastn_mito -- bash -c "tblastn -query '${params.mito_proteins}' -db asm_db -evalue 1e-10 -max_target_seqs 50 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > mito_hits.tsv 2>tblastn.log"
  mv ${meta.id}.nuclear.fasta pre_split.fasta
  ff_run mito_extract -- mito_extract.py --sample "${meta.id}" --assembly pre_split.fasta --tblastn mito_hits.tsv \
      --max-len ${params.mito_max_len} --out-mito ${meta.id}.mito.fasta --out-nuclear ${meta.id}.nuclear.fasta --json mito_extract.json
  python3 - <<'PY'
  import json
  d = json.load(open("${meta.id}.decontam.json")); m = json.load(open("mito_extract.json"))
  d.update(mito_contigs=m["mito_contigs"], mito_bp=m["mito_bp"], mito_candidates=m["candidates"])
  json.dump(d, open("${meta.id}.decontam.json", "w"), indent=2)
  PY
  ff_finalize
  """
  stub:
  // --stub_verdicts 'ID=non_fungal,ID2=human' lets the DAG tests exercise the gate
  def forced = (params.stub_verdicts ?: '').toString().split(',').collect { kv -> kv.split('=') }.findAll { kv -> kv.size() == 2 }.collectEntries { kv -> [(kv[0].trim()): kv[1].trim()] }
  def verdict = forced[meta.id] ?: 'fungal'
  "touch ${meta.id}.nuclear.fasta ${meta.id}.mito.fasta; echo '{\"sample\":\"${meta.id}\",\"stage\":\"decontam\",\"verdict\":\"${verdict}\"}' > ${meta.id}.decontam.json"
}
