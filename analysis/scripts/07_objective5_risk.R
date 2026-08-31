#!/usr/bin/env Rscript
# Objective 5 — drivers of risk. Mixed-effects logistic models (facility random effect)
# of AF-resistance / priority-pathogen carriage across environment strata (mirrors bacterial obj5).
suppressMessages({ library(dplyr); library(readr); library(ggplot2) })
here <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
source(file.path(here, "_theme.R")); D <- ff_dirs()
o <- file.path(D$output, "objective5"); f <- file.path(D$figures, "objective5")
dir.create(o, showWarnings = FALSE, recursive = TRUE); dir.create(f, showWarnings = FALSE, recursive = TRUE)
d <- suppressMessages(read_csv(file.path(D$output, "analysis_dataset.csv"), show_col_types = FALSE)) |>
  filter(!is.na(compartment)) |> mutate(env_class = factor(env_class))
fit <- tryCatch({
  library(lme4)
  m <- glmer(has_af_resistance ~ env_class + (1 | facility), data = d, family = binomial)
  broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
}, error = function(e) {
  # fall back to a simple GLM when n is small / singular
  m <- tryCatch(glm(has_af_resistance ~ env_class, data = d, family = binomial), error = function(e) NULL)
  if (is.null(m)) NULL else broom::tidy(m, conf.int = TRUE, exponentiate = TRUE)
})
if (!is.null(fit)) { readr::write_csv(fit, file.path(o, "risk_model_odds_ratios.csv"))
  message(sprintf("[obj5] risk model fit (%d terms)", nrow(fit)))
} else message("[obj5] insufficient data to fit a risk model yet")
