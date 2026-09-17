# fungiforge — comparative analysis layer (Layer 2)

Objective scripts that turn the per-isolate `master_fungi` rows from Layer 1 into the
air/human/surface **One Health** comparative story, mirroring the bacterial `r-analysis/`.

```bash
# after the Nextflow pipeline has produced results/04_summary/master_fungi.tsv:
analysis/run_downstream.sh results/ [sample_metadata.csv]
```

Inputs from Layer 1: the merged, schema-validated **`04_summary/master_fungi.tsv`** (stage 17;
per-isolate `*.master.tsv` rows are the fallback for older runs) and, when the run-level cohort
stage produced them, **`cohort/clonal_groups.tsv`**, `species_clusters.tsv` and
`bgc_families.tsv`, which `01_merge_metadata.R` copies next to the analysis table so
objective 4 reports genome-level clonal groups and objective 7 the cross-compartment gene-cluster
families instead of candidate pairs.

Flow: `01_merge_metadata.R` (gather the master table + cohort tables) → `02_recover_metadata.R`
(→ `output/analysis_dataset.csv`, the table every objective reads) → objectives:

| Script | Objective |
|---|---|
| `03_objective1_composition.R` | fungal community composition & identification (richness, PERMANOVA, priority pathogens) |
| `04_objective2_resistome.R`  | antifungal resistome by compartment (azole/echinocandin prevalence, cyp51A TR carriage) |
| `05_objective3_mobile.R`     | transposable-element load & mycovirus carriage |
| `06_objective4_transmission.R` | cross-compartment clonality (candidate transmission; SNP/phylo confirm) |
| `07_objective5_risk.R`       | mixed-effects risk drivers (which settings concentrate resistant/priority fungi) |
| `08_objective6_virulence.R`  | virulence & resistance×virulence convergence (azole-R *A. fumigatus* from air) |
| `09_objective7_connectivity.R` | shared BGC / mobile-unit connectivity across compartments (BiG-SCAPE + skani) |
| `10_objective8_novelty.R`    | candidate novel species + novel BGCs / antibiotic-production potential |

Conventions match the bacterial study: self-resolving `../output/objectiveN` + `../figures/objectiveN`,
`FigureN_*` PNG+PDF at 300 dpi, shared `_theme.R` (AIR/HUMAN/SURFACE palette), `set.seed(20260729)`.
Outputs are the substrate for the manuscript + seminar deck (built with the `scideck` engine).
