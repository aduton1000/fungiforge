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
  # and --genemark_key the licence file, which normally sits BESIDE that directory rather than
  # inside it; the profiles bind both. GeneMark reads the key from \$HOME/.gm_key, so HOME is a
  # per-task directory here: the container may provide no HOME, one this UID cannot write, or one
  # shared by every concurrent task.
  # GENEMARK_PATH is exported only once the key is actually in place. funannotate picks GeneMark up
  # from that variable, so exporting it first turned every licensing problem into a failed isolate
  # ("License key '.gm_key' not found") instead of the intended fall back to Augustus alone.
  export HOME="\$PWD/home"; mkdir -p "\$HOME"
  GENEMARK=no; GM_WHY=""
  if [ -z "${gm_dir}" ]; then
    GM_WHY="no --genemark_dir (unpacked GeneMark-ES directory)"
  elif [ ! -x "${gm_dir}/gmes_petap.pl" ]; then
    GM_WHY="${gm_dir}/gmes_petap.pl not found or not executable inside the container"
  elif [ -z "${params.genemark_key ?: ''}" ]; then
    GM_WHY="no --genemark_key (GeneMark-ES will not run without its licence file)"
  elif [ ! -s "${params.genemark_key ?: '/nonexistent'}" ]; then
    GM_WHY="--genemark_key is empty or not visible inside the container: ${params.genemark_key ?: ''}"
  elif ! cp "${params.genemark_key ?: '/nonexistent'}" "\$HOME/.gm_key"; then
    GM_WHY="could not place the licence key at \$HOME/.gm_key"
  else
    GENEMARK=yes
  fi
  if [ "\$GENEMARK" = yes ]; then
    export GENEMARK_PATH="${gm_dir}"; export PATH="${gm_dir}:\$PATH"
    ff_version genemark -- bash -c "gmes_petap.pl 2>&1 | grep -m1 -iE 'version|GeneMark' || echo unknown"
  else
    ff_skip genemark "\$GM_WHY — ab-initio prediction runs without GeneMark"
  fi
  # Species-aware training start (AUGUSTUS_SPECIES / BUSCO_DB), verified against what is staged.
  # The helper prints both assignments; if it cannot run, `eval ""` succeeds silently and the
  # failure would only surface far below as an unbound variable naming neither cause nor helper.
  # So check it here and fall back to the documented defaults instead.
  if TRAIN="\$(annotation_training.py --species ${species} --map '${params.annotation_training_map}' \\
             --augustus-config "\$AUGUSTUS_CONFIG_PATH" --funannotate-db "\$FUNANNOTATE_DB" \\
             --default-species ${params.funannotate_seed} --default-busco ${params.funannotate_busco} \\
             --json training.json)"; then
    eval "\$TRAIN"
  else
    ff_skip annotation_training "annotation_training.py failed; falling back to Augustus ${params.funannotate_seed} / BUSCO set ${params.funannotate_busco}"
  fi
  : "\${AUGUSTUS_SPECIES:=${params.funannotate_seed}}"; : "\${BUSCO_DB:=${params.funannotate_busco}}"
  # funannotate requires short (<=16 char), space-free headers: `funannotate sort` renames them;
  # a description-stripping sed is the fallback.
  ff_run funannotate_sort --optional -- bash -c "funannotate sort -i ${masked} -o clean.fasta -b contig --minlen ${params.min_contig_len} 2>sort.log"
  if [ "\$FF_RC" -ne 0 ]; then sed '/^>/ s/[[:space:]].*//' ${masked} > clean.fasta; fi
  # One definition, used for the GeneMark attempt and for the fallback below.
  run_predict(){ funannotate predict -i clean.fasta -o fun -s "${meta.id}" \\
      --cpus ${task.cpus} --busco_db "\$BUSCO_DB" --busco_seed_species "\$AUGUSTUS_SPECIES" ${params.predict_extra ?: ''}; }
  # A licence that is present but rejected (expired, wrong host) only shows up once prediction is
  # under way, hours in. Retry once without GeneMark rather than losing the isolate; the stage is
  # recorded partial and says why.
  if [ "\$GENEMARK" = yes ]; then
    ff_run funannotate_predict --optional -- run_predict
    if [ "\$FF_RC" -ne 0 ]; then
      ff_skip genemark "funannotate predict failed with GeneMark enabled — retried without it; check the licence key (expiry) and ${gm_dir}/gmes_petap.pl"
      unset GENEMARK_PATH; GENEMARK=no; rm -rf fun
      ff_run funannotate_predict_no_genemark -- run_predict
    fi
  else
    ff_run funannotate_predict -- run_predict
  fi
  cp -r fun/predict_results predict_results
  cp predict_results/*.proteins.fa ${meta.id}.predicted.proteins.faa
  NPROT=\$(grep -c '^>' ${meta.id}.predicted.proteins.faa || true)
  ff_run protein_count -- test "\${NPROT:-0}" -ge 1
  python3 - <<PY
  import json
  import os
  t = json.load(open("training.json")) if os.path.exists("training.json") else {
      "augustus_species": "${params.funannotate_seed}", "busco_db": "${params.funannotate_busco}",
      "basis": "default (training helper did not run)", "notes": []}
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
