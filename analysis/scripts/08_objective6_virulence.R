#!/usr/bin/env Rscript
# Objective 6 — virulence & convergence. Virulence/secretome load and the convergence of
# resistance + virulence (headline analog: azole-resistant A. fumigatus from air).
suppressMessages({ library(dplyr); library(readr); library(ggplot2); library(stringr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective6"); dir.create(o, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  filter(!is.na(compartment))
# resistant AND priority-pathogen (convergence)
conv <- d |> filter(has_af_resistance & priority) |>
  select(sample, species, compartment, facility, resistant_classes, cyp51A_TR)
write_csv(conv, file.path(o, "resistance_virulence_convergence.csv"))
message(sprintf("[obj6] %d resistant priority-pathogen isolates (convergence); e.g. azole-R A. fumigatus from air",
                nrow(conv)))
