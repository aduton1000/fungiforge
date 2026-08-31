#!/usr/bin/env Rscript
# Objective 3 — mobile & repeat elements across compartments: TE load, mycovirus carriage.
suppressMessages({ library(dplyr); library(readr); library(ggplot2) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective3"); f <- file.path(D$figures, "objective3")
dir.create(o, showWarnings = FALSE, recursive = TRUE); dir.create(f, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  mutate(compartment = factor(compartment, levels = COMP_LEVELS)) |> filter(!is.na(compartment)) |>
  mutate(te_percent = suppressWarnings(as.numeric(te_percent)),
         n_mycovirus = suppressWarnings(as.integer(n_mycovirus)))
s <- d |> group_by(compartment) |>
  summarise(mean_te_pct = mean(te_percent, na.rm = TRUE),
            mycovirus_carriers = sum(n_mycovirus > 0, na.rm = TRUE), n = n(), .groups = "drop")
write_csv(s, file.path(o, "mobile_by_compartment.csv"))
p <- ggplot(s, aes(compartment, mean_te_pct, fill = compartment)) + geom_col(width = 0.66) +
  scale_fill_manual(values = COMPARTMENT_COLS, guide = "none") +
  labs(title = "Objective 3 — transposable-element load by compartment", x = NULL, y = "mean TE % of genome")
ff_save(p, "Figure3_mobile", w = 7, h = 5, dir = f)
message(sprintf("[obj3] mean TE%% computed; mycovirus carriers total %d", sum(s$mycovirus_carriers)))
