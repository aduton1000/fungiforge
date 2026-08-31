#!/usr/bin/env Rscript
# Objective 8 — novelty & discovery. Candidate novel species + novel BGCs / antibiotic-
# production potential (the discovery angle: "do any of these fungi make new compounds").
suppressMessages({ library(dplyr); library(readr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective8"); dir.create(o, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE))
novel <- d |> filter(is_novel) |> select(sample, species, compartment, facility, novelty, n_bgc)
write_csv(novel, file.path(o, "candidate_novel_species.csv"))
bgc_rich <- d |> arrange(desc(suppressWarnings(as.numeric(n_bgc)))) |>
  select(sample, species, compartment, n_bgc) |> head(20)
write_csv(bgc_rich, file.path(o, "top_bgc_producers.csv"))
message(sprintf("[obj8] %d candidate novel-species isolates; top BGC producers tabulated (novel-cluster call via BiG-SCAPE vs MIBiG)",
                nrow(novel)))
