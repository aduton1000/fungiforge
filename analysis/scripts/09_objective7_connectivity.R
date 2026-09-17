#!/usr/bin/env Rscript
# Objective 7 — connectivity. Shared near-identical BGC / mobile units across compartments
# and species (fungal analog of the bacterial mobile-gene connectivity; BiG-SCAPE families
# + skani units feed this). Scaffold: summarise shared BGC types across compartments.
suppressMessages({ library(dplyr); library(readr); library(tidyr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective7"); dir.create(o, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  filter(!is.na(compartment))
s <- d |> group_by(compartment) |> summarise(mean_bgc = mean(suppressWarnings(as.numeric(n_bgc)), na.rm = TRUE),
                                             n = n(), .groups = "drop")
write_csv(s, file.path(o, "bgc_by_compartment.csv"))

# W3.2/W3.3: gene-cluster families from the cohort stage (BiG-SCAPE, else shared-protein clustering)
fam_file <- file.path(D$output, "bgc_families.tsv")
if (file.exists(fam_file)) {
  fam <- suppressMessages(read_tsv(fam_file, show_col_types = FALSE)) |>
    left_join(d |> select(sample, species, compartment, facility), by = "sample")
  shared <- fam |> group_by(family) |>
    summarise(n_regions = n(), n_isolates = n_distinct(sample),
              species = paste(sort(unique(na.omit(species))), collapse = ";"), n_species = n_distinct(na.omit(species)),
              compartments = paste(sort(unique(as.character(na.omit(compartment)))), collapse = ";"),
              n_compartments = n_distinct(na.omit(compartment)),
              products = paste(sort(unique(na.omit(products))), collapse = ";"),
              mibig = paste(sort(unique(na.omit(mibig_best[mibig_best != ""]))), collapse = ";"), .groups = "drop") |>
    filter(n_isolates > 1) |> arrange(desc(n_compartments), desc(n_isolates))
  write_csv(shared, file.path(o, "shared_bgc_families.csv"))
  cross <- shared |> filter(n_compartments > 1)
  write_csv(cross, file.path(o, "cross_compartment_bgc_families.csv"))
  message(sprintf("[obj7] %d gene-cluster families shared by >1 isolate; %d span >1 compartment; %d cross species",
                  nrow(shared), nrow(cross), sum(shared$n_species > 1)))
} else {
  message("[obj7] no cohort/bgc_families.tsv — per-compartment BGC counts only (enable the cohort stage)")
}
