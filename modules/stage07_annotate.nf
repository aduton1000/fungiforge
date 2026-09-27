// Stage 07 — functional annotation (Funannotate annotate; W2.5 split). Takes the prediction of
// stage 07a plus the eggNOG-mapper table (07b) and the InterProScan XML (07c) when present, and
// runs funannotate annotate (Pfam, dbCAN, MEROPS, BUSCO, UniProt/SwissProt, plus --eggnog /
// --iprscan when given). Emits proteins + GBK (consumed by stages 09-13) and annotate.json with
// the annotation coverage (bin/annotate_stats.py).
process ANNOTATE {
  tag { meta.id }
  label 'annotate'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy'
  input:  tuple val(meta), path(predict_results), path(eggnog), path(iprscan)
  output: tuple val(meta), path("${meta.id}.proteins.faa"), emit: proteins
          tuple val(meta), path("${meta.id}.gbk"),          emit: gbk
          tuple val(meta), path("${meta.id}.annotate.json"),emit: json
  script:
  """
  source ff_status.sh; ff_init annotate "${meta.id}" ${meta.id}.annotate.json
  ff_version funannotate -- funannotate --version
  ff_version diamond -- diamond version
  export FUNANNOTATE_DB="${params.funannotate_db ?: params.data_dir + '/funannotate'}"
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  # funannotate annotate reads <dir>/predict_results and writes annotate_misc / annotate_results:
  # a private copy so the staged prediction (a symlink into stage 07a's task) is never written to.
  mkdir -p fun && cp -r ${predict_results}/. fun/predict_results/
  EGG=""; if [ -s ${eggnog} ]; then EGG="--eggnog ${eggnog}"; EGGS=yes; else ff_skip eggnog "no eggNOG-mapper table (stage 07b skipped or off)"; EGGS=no; fi
  IPS=""; if [ -s ${iprscan} ]; then IPS="--iprscan ${iprscan}"; IPSS=yes; else ff_skip interproscan "no InterProScan XML (stage 07c skipped or off)"; IPSS=no; fi
  ff_run funannotate_annotate -- funannotate annotate -i fun --cpus ${task.cpus} \$EGG \$IPS ${params.annotate_extra ?: ''}
  if ls fun/annotate_results/*.proteins.fa >/dev/null 2>&1; then cp fun/annotate_results/*.proteins.fa ${meta.id}.proteins.faa
  else cp fun/predict_results/*.proteins.fa ${meta.id}.proteins.faa; fi
  if ls fun/annotate_results/*.gbk >/dev/null 2>&1; then cp fun/annotate_results/*.gbk ${meta.id}.gbk
  else cp fun/predict_results/*.gbk ${meta.id}.gbk; fi
  # A 0-protein result would silently poison resistance/BGC downstream — treat it as a failure.
  NPROT=\$(grep -c '^>' ${meta.id}.proteins.faa || true)
  ff_run protein_count -- test "\${NPROT:-0}" -ge 1
  ANN=\$(ls fun/annotate_results/*.annotations.txt 2>/dev/null | head -1)
  GM=\$(grep -o '"genemark": *"[a-z]*"' ${predict_results}/../*.predict.json 2>/dev/null | head -1 | grep -o '[a-z]*"\$' | tr -d '"')
  ff_run annotate_stats -- annotate_stats.py --sample "${meta.id}" --annotations "\${ANN:-none}" --proteins ${meta.id}.proteins.faa \\
      --eggnog \$EGGS --interproscan \$IPSS --genemark "\${GM:-unknown}" --out ${meta.id}.annotate.json
  ff_finalize
  """
  stub:
  "touch ${meta.id}.proteins.faa ${meta.id}.gbk; echo '{\"sample\":\"${meta.id}\",\"stage\":\"annotate\"}' > ${meta.id}.annotate.json"
}
