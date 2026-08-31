#!/usr/bin/env Rscript
# Objective 2 — antifungal resistome across compartments. Prevalence of AF resistance,
# by drug class (azole/echinocandin), cyp51A TR carriage, known-mutation load; Fisher tests.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2); library(stringr) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o2 <- file.path(D$output, "objective2"); dir.create(o2, showWarnings = FALSE, recursive = TRUE)
f2 <- file.path(D$figures, "objective2"); dir.create(f2, showWarnings = FALSE, recursive = TRUE)

d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  mutate(compartment = factor(compartment, levels = COMP_LEVELS)) |> filter(!is.na(compartment))

wilson <- function(x, n) if (n == 0) c(NA, NA, NA) else {
  p <- x/n; z <- 1.96; den <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/den; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(p, max(0, ctr-hw), min(1, ctr+hw))
}

# resistance prevalence by compartment
prev <- d |> group_by(compartment) |>
  summarise(n = n(), any_R = sum(has_af_resistance, na.rm = TRUE),
            azole_R = sum(azole_R, na.rm = TRUE), echino_R = sum(echino_R, na.rm = TRUE),
            multi_R = sum(multi_af_R, na.rm = TRUE), .groups = "drop") |>
  rowwise() |> mutate(p_any = wilson(any_R, n)[1], lo = wilson(any_R, n)[2], hi = wilson(any_R, n)[3]) |> ungroup()
write_csv(prev, file.path(o2, "resistance_prevalence_by_compartment.csv"))

fish <- tryCatch(fisher.test(table(d$compartment, d$has_af_resistance))$p.value, error = function(e) NA)

# cyp51A TR carriage (A. fumigatus azole mechanism)
tr <- d |> filter(!is.na(cyp51A_TR)) |> count(compartment, cyp51A_TR)
write_csv(tr, file.path(o2, "cyp51A_TR_carriage.csv"))

pl <- prev |> select(compartment, azole = azole_R, echinocandin = echino_R, n) |>
  pivot_longer(c(azole, echinocandin), names_to = "class", values_to = "count") |>
  mutate(frac = count / n)
p <- ggplot(pl, aes(compartment, frac, fill = class)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = count), position = position_dodge(0.75), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(azole = "#8C3B4A", echinocandin = "#1B5E20"), name = "drug class") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Objective 2 — antifungal resistance by compartment",
       subtitle = sprintf("any-resistance Fisher p = %s",
                          ifelse(is.na(fish), "NA", format.pval(fish, digits = 2))),
       x = NULL, y = "% of isolates")
ff_save(p, "Figure2_resistome", w = 7.5, h = 5, dir = f2)
message(sprintf("[obj2] AF-resistant: %d/%d; azole %d; echinocandin %d; cyp51A-TR carriers %d",
                sum(d$has_af_resistance, na.rm = TRUE), nrow(d),
                sum(d$azole_R, na.rm = TRUE), sum(d$echino_R, na.rm = TRUE), nrow(tr)))
