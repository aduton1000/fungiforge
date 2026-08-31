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
  export FUNANNOTATE_DB="${params.funannotate_db ?: params.data_dir + '/funannotate'}"
  # Augustus (funannotate's BUSCO self-training) must WRITE species params into
  # AUGUSTUS_CONFIG_PATH; the container default (/usr/share/augustus/config) is root-owned
  # (read-only to our host UID) so Augustus fails silently -> 0 BUSCOs -> predict aborts.
  # funannotate derives AUGUSTUS_BASE = dirname(CONFIG) ONLY when basename(CONFIG)=='config',
  # then needs BASE/scripts/*.pl and BASE/bin/bam2hints. So mirror the layout: a writable
  # 'config' copy with scripts/ and bin/bam2hints symlinked to the real (read-only) ones.
  export AUGUSTUS_CONFIG_PATH=/tmp/augustus/config
  if [ ! -d "\$AUGUSTUS_CONFIG_PATH/species" ]; then
    mkdir -p /tmp/augustus/bin
    cp -r /usr/share/augustus/config /tmp/augustus/config 2>/dev/null || true
    ln -sf /usr/share/augustus/scripts /tmp/augustus/scripts 2>/dev/null || true
    ln -sf "\$(command -v bam2hints || echo /usr/bin/bam2hints)" /tmp/augustus/bin/bam2hints 2>/dev/null || true
  fi
  ${ params.genemark_key ? "cp ${params.genemark_key} ~/.gm_key || true" : "" }
  # funannotate requires short (<=16 char), space-free FASTA headers. Polypolish appends
  # " polypolish" to every contig (>contig_30 polypolish) which blows the 16-char limit and
  # aborts predict. funannotate 'sort' renames headers (contig_1..N, longest first) + is the
  # documented preprocessing step; sed-strip of the description is the fallback.
  funannotate sort -i ${masked} -o clean.fasta -b contig 2>sort.log || sed '/^>/ s/[[:space:]].*//' ${masked} > clean.fasta
  funannotate predict -i clean.fasta -o fun -s "${meta.id}" \\
      --cpus ${task.cpus} --busco_db ${params.funannotate_busco} || true
  EGG=""; [ -n "\$(ls -A ${params.data_dir}/eggnog 2>/dev/null)" ] && EGG="--eggnog ${params.data_dir}/eggnog"
  funannotate annotate -i fun --cpus ${task.cpus} ${ips} \$EGG || true
  cp fun/predict_results/*.proteins.fa ${meta.id}.proteins.faa 2>/dev/null || \\
     cp fun/annotate_results/*.proteins.fa ${meta.id}.proteins.faa 2>/dev/null || touch ${meta.id}.proteins.faa
  cp fun/annotate_results/*.gbk ${meta.id}.gbk 2>/dev/null || cp fun/predict_results/*.gbk ${meta.id}.gbk 2>/dev/null || touch ${meta.id}.gbk
  # Fail LOUD on an empty annotation: a 0-protein result silently poisons resistance/BGC
  # downstream (they find nothing and report false 'not_detected'). grep -c prints "0" and
  # exits 1 on no match, so guard the substitution, don't append a second 0.
  NPROT=\$(grep -c '^>' ${meta.id}.proteins.faa 2>/dev/null || true); NPROT=\${NPROT:-0}
  if [ "\$NPROT" -lt 1 ]; then
    echo "ERROR: funannotate produced 0 proteins for ${meta.id} — annotation failed (see fun/logfiles + sort.log)." >&2
    exit 1
  fi
  printf '{"sample":"%s","stage":"annotate","proteins":"%s.proteins.faa","n_proteins":%s,"interproscan":%s}\\n' \\
    "${meta.id}" "${meta.id}" "\$NPROT" "${params.run_interproscan.toString()}" > ${meta.id}.annotate.json
  """
  stub:
  "touch ${meta.id}.proteins.faa ${meta.id}.gbk ${meta.id}.annotate.json"
}
