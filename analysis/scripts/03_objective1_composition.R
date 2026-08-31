#!/usr/bin/env Rscript
# Objective 1 — fungal community composition & identification across compartments.
# Richness, compartment PERMANOVA (species presence/absence), priority-pathogen carriage.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2); library(stringr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o1 <- file.path(D$output, "objective1"); dir.create(o1, showWarnings = FALSE, recursive = TRUE)
f1 <- file.path(D$figures, "objective1"); dir.create(f1, showWarnings = FALSE, recursive = TRUE)

d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  mutate(compartment = factor(compartment, levels = COMP_LEVELS)) |> filter(!is.na(compartment))

# richness per compartment
rich <- d |> group_by(compartment) |>
  summarise(n_isolates = n(), n_species = n_distinct(species), .groups = "drop")
write_csv(rich, file.path(o1, "richness_by_compartment.csv"))

# compartment PERMANOVA on species presence/absence (needs vegan + >1 species/compartment)
perm <- tryCatch({
  library(vegan)
  mat <- d |> mutate(present = 1L) |> distinct(sample, species, compartment, present) |>
    pivot_wider(names_from = species, values_from = present, values_fill = 0)
  X <- as.matrix(mat |> select(-sample, -compartment))
  if (nrow(X) >= 4 && ncol(X) >= 2) {
    a <- adonis2(X ~ compartment, data = mat, method = "jaccard", permutations = 999)
    tibble(term = "compartment", R2 = a$R2[1], p = a$`Pr(>F)`[1])
  } else tibble(term = "compartment", R2 = NA, p = NA)
}, error = function(e) tibble(term = "compartment", R2 = NA, p = NA, note = conditionMessage(e)))
write_csv(perm, file.path(o1, "permanova_community.csv"))

# priority fungal pathogens by compartment
prio <- d |> filter(priority) |> count(species, compartment) |> arrange(desc(n))
write_csv(prio, file.path(o1, "priority_pathogens.csv"))

p <- ggplot(rich, aes(compartment, n_species, fill = compartment)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%d spp\n(%d iso)", n_species, n_isolates)), vjust = -0.2, size = 3) +
  scale_fill_manual(values = COMPARTMENT_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Objective 1 — fungal richness by compartment",
       subtitle = sprintf("PERMANOVA R2 = %s, p = %s",
                          ifelse(is.na(perm$R2[1]), "NA", sprintf("%.3f", perm$R2[1])),
                          ifelse(is.na(perm$p[1]),  "NA", sprintf("%.3f", perm$p[1]))),
       x = NULL, y = "distinct species")
ff_save(p, "Figure1_composition", w = 7, h = 5, dir = f1)
message(sprintf("[obj1] %d species across %d isolates; %d priority-pathogen isolates",
                sum(rich$n_species), sum(rich$n_isolates), sum(d$priority, na.rm = TRUE)))
