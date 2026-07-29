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
# placeholder connectivity summary from per-isolate BGC counts; the real edge list comes
# from BiG-SCAPE GCF membership + skani mobile-unit ANI (built in a bigscape aux step).
s <- d |> group_by(compartment) |> summarise(mean_bgc = mean(suppressWarnings(as.numeric(n_bgc)), na.rm = TRUE),
                                             n = n(), .groups = "drop")
write_csv(s, file.path(o, "bgc_by_compartment.csv"))
message("[obj7] BGC-per-compartment summarised; run BiG-SCAPE aux for the cross-compartment GCF network")
