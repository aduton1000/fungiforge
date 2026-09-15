// Stage 07 — eukaryotic gene prediction + functional annotation (Funannotate).
// predict (GeneMark/Augustus/SNAP/GlimmerHMM via EVM) -> annotate (eggNOG, Pfam,
// dbCAN, MEROPS, antiSMASH handoff; InterProScan optional). Emits proteins + GBK.
process ANNOTATE {
  tag { meta.id }
  label 'annotate'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy'
  input:  tuple val(meta), path(masked)
  output: tuple val(meta), path("${meta.id}.proteins.faa"), emit: proteins
          tuple val(meta), path("${meta.id}.gbk"),          emit: gbk
          tuple val(meta), path("${meta.id}.annotate.json"),emit: json
  script:
  // funannotate 'annotate' runs InterProScan only when given --iprscan (an XML from a
  // separate run); there is NO --no-interpro flag (passing it aborts annotate). So we
  // simply omit iprscan unless a pre-computed XML is supplied via --iprscan_xml.
  def ips = (params.run_interproscan && params.iprscan_xml) ? "--iprscan ${params.iprscan_xml}" : ''
  // --busco_db takes a funannotate BUSCO set name (dikarya|fungi|…), NOT a BUSCO lineage
  // like 'fungi_odb10' and NOT an Augustus species. 'dikarya' ships with the funannotate
  // DB and covers Asco/Basidiomycota. --busco_seed_species (Augustus species) is left at
  // funannotate's default (anidulans) for organism-agnostic self-training.
  """
  source ff_status.sh; ff_init annotate "${meta.id}" ${meta.id}.annotate.json
  ff_version funannotate -- funannotate --version
  ff_version augustus -- augustus --version
  ff_version genemark -- gmes_petap.pl
  ff_version emapper -- emapper.py --version
  ff_version diamond -- diamond version
  export FUNANNOTATE_DB="${params.funannotate_db ?: params.data_dir + '/funannotate'}"
  # Scratch stays in the task directory: the container's /tmp is the node's tmpfs (RAM-backed,
  # shared by every task on the node). funannotate writes the protein-to-genome alignments
  # (exonerate, tens of GB for a 500k-protein database) under TMPDIR, and a batch of concurrent
  # annotate tasks filled it ("No space left on device"). The work dir is on the shared disk.
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  # Augustus (funannotate's BUSCO self-training) must WRITE species params into
  # AUGUSTUS_CONFIG_PATH; the container default (/usr/share/augustus/config) is root-owned
  # (read-only to our host UID) so Augustus fails silently -> 0 BUSCOs -> predict aborts.
  # funannotate derives AUGUSTUS_BASE = dirname(CONFIG) ONLY when basename(CONFIG)=='config',
  # then needs BASE/scripts/*.pl and BASE/bin/bam2hints. So mirror the layout: a writable
  # 'config' copy with scripts/ and bin/bam2hints symlinked to the real (read-only) ones.
  # (per-task copy under the work dir, not the node-shared /tmp: concurrent tasks must not race)
  export AUGUSTUS_CONFIG_PATH="\$PWD/augustus/config"
  if [ ! -d "\$AUGUSTUS_CONFIG_PATH/species" ]; then
    mkdir -p "\$PWD/augustus/bin"
    [ -d /usr/share/augustus/config ]  && cp -r /usr/share/augustus/config "\$PWD/augustus/config"
    [ -d /usr/share/augustus/scripts ] && ln -sf /usr/share/augustus/scripts "\$PWD/augustus/scripts"
    B2H=\$(command -v bam2hints || true); [ -n "\$B2H" ] && ln -sf "\$B2H" "\$PWD/augustus/bin/bam2hints"
  fi
  ${ params.genemark_key ? "cp ${params.genemark_key} ~/.gm_key" : "ff_skip genemark 'no --genemark_key given (GeneMark-ES ab-initio prediction skipped by funannotate)'" }
  # funannotate requires short (<=16 char), space-free FASTA headers. Polypolish appends
  # " polypolish" to every contig (>contig_30 polypolish) which blows the 16-char limit and
  # aborts predict. funannotate 'sort' renames headers (contig_1..N, longest first) + is the
  # documented preprocessing step; sed-strip of the description is the fallback (partial).
  ff_run funannotate_sort --optional -- bash -c "funannotate sort -i ${masked} -o clean.fasta -b contig --minlen ${params.min_contig_len} 2>sort.log"
  if [ "\$FF_RC" -ne 0 ]; then sed '/^>/ s/[[:space:]].*//' ${masked} > clean.fasta; fi
  ff_run funannotate_predict -- funannotate predict -i clean.fasta -o fun -s "${meta.id}" \\
      --cpus ${task.cpus} --busco_db ${params.funannotate_busco}
  EGG=""; if [ -n "\$(ls -A ${params.data_dir}/eggnog 2>/dev/null)" ]; then EGG="--eggnog ${params.data_dir}/eggnog"; else ff_skip eggnog "eggNOG database not staged under --data_dir"; fi
  ff_run funannotate_annotate -- funannotate annotate -i fun --cpus ${task.cpus} ${ips} \$EGG
  if ls fun/annotate_results/*.proteins.fa >/dev/null 2>&1; then cp fun/annotate_results/*.proteins.fa ${meta.id}.proteins.faa
  else cp fun/predict_results/*.proteins.fa ${meta.id}.proteins.faa; fi
  if ls fun/annotate_results/*.gbk >/dev/null 2>&1; then cp fun/annotate_results/*.gbk ${meta.id}.gbk
  else cp fun/predict_results/*.gbk ${meta.id}.gbk; fi
  # A 0-protein result would silently poison resistance/BGC downstream — treat it as a failure.
  NPROT=\$(grep -c '^>' ${meta.id}.proteins.faa || true)
  ff_run protein_count -- test "\${NPROT:-0}" -ge 1
  printf '{"sample":"%s","stage":"annotate","proteins":"%s.proteins.faa","n_proteins":%s,"interproscan":%s}\\n' \\
    "${meta.id}" "${meta.id}" "\$NPROT" "${(params.run_interproscan && params.iprscan_xml) ? 'true' : 'false'}" > ${meta.id}.annotate.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.proteins.faa ${meta.id}.gbk; echo '{\"sample\":\"${meta.id}\",\"stage\":\"annotate\"}' > ${meta.id}.annotate.json"
}
