// Stage 16 (run-level) — cohort phylogenomics and clonality (W3.1). Runs once after every isolate
// that passed the gates: skani all-vs-all ANI -> species clusters; single-copy BUSCO orthologs
// shared by the cohort -> MAFFT -> gap-trimmed supermatrix -> IQ-TREE 3 (or FastTree);
// within each species cluster, assembly-based SNP distances to the most complete member
// (minimap2 asm5 + paftools.js call, SNVs inside the core aligned in every member) and clonal
// groups at --clonal_snp_threshold; gene-cluster families across the isolates' antiSMASH regions
// (W3.2, bin/bgc_families.py: BiG-SCAPE 2 with the funannotate database's Pfam when available,
// else DIAMOND shared-protein clustering). Outputs go to <outdir>/cohort/. --skip_cohort removes it.
process COHORT {
  tag 'run'
  label 'cohort'
  publishDir { "${params.outdir}" }, mode: 'copy', pattern: "cohort/*"   // -> <outdir>/cohort/
  input:  path(assemblies, stageAs: 'asm/*')
          path(busco,      stageAs: 'busco/*')
          path(species,    stageAs: 'species/*')
          path(qc,         stageAs: 'qc/*')
          path(regions,    stageAs: 'regions/*')
          path(bgc,        stageAs: 'bgc/*')
  output: path('cohort/cohort.json'), emit: json
          path('cohort/*'),           emit: files
  script:
  """
  source ff_status.sh; ff_init cohort run cohort.json --best-effort
  ff_version skani -- skani --version
  ff_version mafft -- bash -c "mafft --version 2>&1 | head -1"
  ff_version iqtree3 -- bash -c "iqtree3 --version 2>&1 | head -1"
  ff_version minimap2 -- minimap2 --version
  ff_run cohort_phylo -- cohort_phylo.py --outdir cohort --threads ${task.cpus} \\
      --assemblies asm/* --busco busco/* --species species/* --qc qc/* \\
      --min-frac ${params.cohort_min_frac} --min-isolates ${params.cohort_min_isolates} --tree ${params.cohort_tree_method} \\
      --snp-threshold ${params.clonal_snp_threshold} --ani-cluster ${params.cohort_ani_cluster} --max-genes ${params.cohort_max_genes}
  # W3.2: gene-cluster families across the cohort from the antiSMASH region GenBanks
  ff_version diamond -- diamond version
  ff_version bigscape -- bash -c "bigscape --version 2>&1 | head -1"
  ff_run bgc_families -- bgc_families.py --outdir cohort --threads ${task.cpus} --regions regions/* --bgc bgc/* \\
      --method ${params.gcf_method} --gcf-cutoff ${params.gcf_cutoff} --data-dir "${params.data_dir ?: ''}" \\
      ${params.funannotate_db ? "--pfam-path ${params.funannotate_db}/Pfam-A.hmm" : ''} \\
      --min-sim ${params.gcf_min_similarity} --min-pident ${params.gcf_min_pident}
  # The stage status rides inside cohort/cohort.json: a run-level stage has no per-sample stage
  # JSON. ff_finalize must run FIRST — it is what writes the status document this block reads.
  ff_finalize
  python3 - <<'PY'
  import json, os
  d = json.load(open("cohort/cohort.json")) if os.path.exists("cohort/cohort.json") else {"stage": "cohort"}
  if os.path.exists("cohort.json"):
      st = json.load(open("cohort.json"))
      for k in ("status", "tools", "skipped_tools", "versions"):
          if k in st:
              d[k] = st[k]
      d["stage_status"] = st.get("status")
  try:
      d["bgc_families"] = {k: v for k, v in json.load(open("cohort/bgc_families.json")).items() if k != "families"}
  except Exception:
      pass
  os.makedirs("cohort", exist_ok=True)
  json.dump(d, open("cohort/cohort.json", "w"), indent=2)
  PY
  """
  stub:
  "mkdir -p cohort; echo '{\"stage\":\"cohort\",\"n_isolates\":0}' > cohort/cohort.json; touch cohort/species_clusters.tsv cohort/bgc_families.tsv"
}
