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
message(sprintf("[obj4] %d species span >1 compartment (candidate transmission; confirm by SNP/phylogeny)",
                nrow(span)))
