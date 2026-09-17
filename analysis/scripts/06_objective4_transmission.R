#!/usr/bin/env Rscript
# Objective 4 — clonal transmission across compartments. Same-species isolates in >1
# compartment are candidate transmission pairs; SNP/phylogenomic confirmation (skani/ska2
# + iqtree) is added when >=2 same-species genomes exist (mirrors bacterial objective 4).
suppressMessages({ library(dplyr); library(readr); library(tidyr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective4"); dir.create(o, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  filter(!is.na(compartment))
# species spanning >1 compartment (candidate cross-compartment transmission)
span <- d |> group_by(species) |>
  summarise(n = n(), compartments = paste(sort(unique(as.character(compartment))), collapse = ";"),
            n_comp = n_distinct(compartment), .groups = "drop") |>
  filter(n_comp > 1) |> arrange(desc(n_comp), desc(n))
write_csv(span, file.path(o, "cross_compartment_species.csv"))

# W3.1/W3.3: the pipeline's cohort stage gives genome-level confirmation. clonal_groups.tsv
# (single linkage at the run's SNP threshold) turns candidates into evidenced transmission links.
cg_file <- file.path(D$output, "clonal_groups.tsv")
if (file.exists(cg_file)) {
  cg <- suppressMessages(read_tsv(cg_file, show_col_types = FALSE))
  linked <- cg |> group_by(cluster, clonal_group) |> filter(n() > 1) |> ungroup() |>
    left_join(d |> select(sample, species, compartment, facility), by = "sample")
  write_csv(linked, file.path(o, "clonal_group_members.csv"))
  cross <- linked |> group_by(cluster, clonal_group) |>
    summarise(n_isolates = n(), species = paste(sort(unique(na.omit(species))), collapse = ";"),
              compartments = paste(sort(unique(as.character(na.omit(compartment)))), collapse = ";"),
              n_compartments = n_distinct(na.omit(compartment)),
              facilities = paste(sort(unique(na.omit(facility))), collapse = ";"),
              n_facilities = n_distinct(na.omit(facility)), .groups = "drop") |>
    arrange(desc(n_compartments), desc(n_isolates))
  write_csv(cross, file.path(o, "clonal_groups_cross_compartment.csv"))
  message(sprintf("[obj4] %d clonal group(s) of >1 isolate; %d span >1 compartment (genome-level, from cohort/clonal_groups.tsv)",
                  nrow(cross), sum(cross$n_compartments > 1)))
} else {
  message("[obj4] no cohort/clonal_groups.tsv — candidate pairs only (run the pipeline with the cohort stage enabled)")
}
message(sprintf("[obj4] %d species span >1 compartment (candidate transmission)", nrow(span)))
