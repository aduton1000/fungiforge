// Stage 07a — gene prediction (Funannotate predict; W2.5). GeneMark-ES joins Augustus, SNAP and
// GlimmerHMM in EVM when the licensed GeneMark directory is supplied (--genemark_dir, bound
// into the container, plus --genemark_key). Training is species-aware: the stage-08 call picks
// the Augustus pre-trained species to seed BUSCO training and the funannotate BUSCO set
// (fungiforge/resources/annotation_training.tsv; bin/annotation_training.py falls back to
// anidulans / dikarya when a mapped choice is not staged). Emits predict_results/ for stage 07
// (funannotate annotate) and the predicted proteins for the eggNOG / InterProScan stages.
process PREDICT {
  tag { meta.id }
  label 'annotate'
  publishDir { "${params.outdir}/${meta.id}/07_annotate" }, mode: 'copy', pattern: "*.json"
  input:  tuple val(meta), path(masked), path(species)
  output: tuple val(meta), path("predict_results"),                 emit: results
          tuple val(meta), path("${meta.id}.predicted.proteins.faa"), emit: proteins
          tuple val(meta), path("${meta.id}.predict.json"),         emit: json
  script:
  def gm_dir = params.genemark_dir ?: ''
  """
  source ff_status.sh; ff_init predict "${meta.id}" ${meta.id}.predict.json
  ff_version funannotate -- funannotate --version
  ff_version augustus -- augustus --version
  export FUNANNOTATE_DB="${params.funannotate_db ?: params.data_dir + '/funannotate'}"
  # Scratch stays in the task directory (the container's /tmp is the node's tmpfs; see stage 07).
  export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"
  # Augustus must WRITE species params into AUGUSTUS_CONFIG_PATH; the container default is
  # read-only to the host UID, so mirror the layout into a writable per-task copy (see stage 07 notes).
  export AUGUSTUS_CONFIG_PATH="\$PWD/augustus/config"
  if [ ! -d "\$AUGUSTUS_CONFIG_PATH/species" ]; then
    mkdir -p "\$PWD/augustus/bin"
    [ -d /usr/share/augustus/config ]  && cp -r /usr/share/augustus/config "\$PWD/augustus/config"
    [ -d /usr/share/augustus/scripts ] && ln -sf /usr/share/augustus/scripts "\$PWD/augustus/scripts"
    B2H=\$(command -v bam2hints || true); [ -n "\$B2H" ] && ln -sf "\$B2H" "\$PWD/augustus/bin/bam2hints"
  fi
  # GeneMark-ES: licensed, not in the image. --genemark_dir is the unpacked gmes_linux_64 directory
  # (bound into the container by the profile); the key goes to \$HOME/.gm_key (HOME is per task).
  GENEMARK=no
  if [ -n "${gm_dir}" ] && [ -x "${gm_dir}/gmes_petap.pl" ]; then
    export GENEMARK_PATH="${gm_dir}"; export PATH="${gm_dir}:\$PATH"
    if [ -n "${params.genemark_key ?: ''}" ] && [ -s "${params.genemark_key ?: ''}" ]; then
      export HOME="\${HOME:-\$PWD}"; mkdir -p "\$HOME" 2>/dev/null || export HOME="\$PWD"
      cp "${params.genemark_key}" "\$HOME/.gm_key" && GENEMARK=yes
      ff_version genemark -- bash -c "gmes_petap.pl 2>&1 | grep -m1 -iE 'version|GeneMark' || echo unknown"
    else
      ff_skip genemark "no --genemark_key: GeneMark-ES needs the licence key next to the binaries"
    fi
  else
    ff_skip genemark "no --genemark_dir (unpacked GeneMark-ES directory): ab-initio prediction runs without GeneMark"
  fi
  # Species-aware training start (AUGUSTUS_SPECIES / BUSCO_DB), verified against what is staged.
  eval "\$(annotation_training.py --species ${species} --map '${params.annotation_training_map}' \\
             --augustus-config "\$AUGUSTUS_CONFIG_PATH" --funannotate-db "\$FUNANNOTATE_DB" \\
             --default-busco ${params.funannotate_busco} --json training.json)"
  # funannotate requires short (<=16 char), space-free headers: `funannotate sort` renames them;
  # a description-stripping sed is the fallback.
  ff_run funannotate_sort --optional -- bash -c "funannotate sort -i ${masked} -o clean.fasta -b contig --minlen ${params.min_contig_len} 2>sort.log"
  if [ "\$FF_RC" -ne 0 ]; then sed '/^>/ s/[[:space:]].*//' ${masked} > clean.fasta; fi
  ff_run funannotate_predict -- funannotate predict -i clean.fasta -o fun -s "${meta.id}" \\
      --cpus ${task.cpus} --busco_db "\$BUSCO_DB" --busco_seed_species "\$AUGUSTUS_SPECIES" ${params.predict_extra ?: ''}
  cp -r fun/predict_results predict_results
  cp predict_results/*.proteins.fa ${meta.id}.predicted.proteins.faa
  NPROT=\$(grep -c '^>' ${meta.id}.predicted.proteins.faa || true)
  ff_run protein_count -- test "\${NPROT:-0}" -ge 1
  python3 - <<PY
  import json
  t = json.load(open("training.json"))
  json.dump({"sample": "${meta.id}", "stage": "predict", "n_proteins": int("\${NPROT:-0}"), "genemark": "\$GENEMARK",
             "augustus_species": t["augustus_species"], "busco_db": t["busco_db"], "training_basis": t["basis"], "training_notes": t["notes"]},
            open("${meta.id}.predict.json", "w"), indent=2)
  PY
  ff_finalize
  """
  stub:
  """
  mkdir -p predict_results; touch predict_results/${meta.id}.gbk ${meta.id}.predicted.proteins.faa
  echo '{"sample":"${meta.id}","stage":"predict","genemark":"no"}' > ${meta.id}.predict.json
  """
}
