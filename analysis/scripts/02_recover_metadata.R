#!/usr/bin/env Rscript
# 02 — standardise into the analysis-ready analysis_dataset.csv that every objective
# reads. Normalises compartment/env, derives resistance/novelty/BGC logicals.
suppressMessages({ library(dplyr); library(readr); library(stringr); library(tidyr) })
source(file.path(dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])), "_theme.R"))
D <- ff_dirs()
m <- suppressMessages(read_tsv(file.path(D$output, "merged_master.tsv"), show_col_types = FALSE))

norm_comp <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- dplyr::recode(x, DOOR = "SURFACE", SWAB = "HUMAN", HAND = "HUMAN", .default = x)
  factor(ifelse(x %in% COMP_LEVELS, x, NA), levels = COMP_LEVELS)
}
env_of <- function(fac) {
  f <- tolower(as.character(fac))
  dplyr::case_when(str_detect(f, "abattoir|farm|food|market") ~ "Food/Agri",
                   str_detect(f, "clinic|hospital|ward|dental") ~ "Clinical",
                   TRUE ~ "Community")
}

d <- m |>
  mutate(compartment = norm_comp(compartment),
         env_class = env_of(facility),
         season = str_to_title(tolower(trimws(as.character(season)))),
         # resistance
         resistant_classes = ifelse(is.na(resistant_classes), "none", resistant_classes),
         n_af_classes = ifelse(resistant_classes == "none", 0L,
                               str_count(resistant_classes, ";") + 1L),
         has_af_resistance = n_af_classes >= 1,
         azole_R  = str_detect(resistant_classes, "azole"),
         echino_R = str_detect(resistant_classes, "echinocandin"),
         multi_af_R = n_af_classes >= 2,
         cyp51A_TR = ifelse(is.na(cyp51A_TR) | cyp51A_TR %in% c("none","NA"), NA, cyp51A_TR),
         # discovery
         is_novel = str_detect(tolower(coalesce(novelty, "")), "novel"),
         n_bgc = suppressWarnings(as.integer(n_bgc)),
         n_known_af_mutations = suppressWarnings(as.integer(n_known_af_mutations)),
         priority = str_detect(species, regex("Candida auris|Aspergillus fumigatus|Cryptococcus|Candida albicans|Candida glabrata|Nakaseomyces|Fusarium", ignore_case = TRUE)))

write_csv(d, file.path(D$output, "analysis_dataset.csv"))
message(sprintf("[02] analysis_dataset.csv: %d isolates | %d compartment-known | %d AF-resistant | %d priority | %d novel",
                nrow(d), sum(!is.na(d$compartment)), sum(d$has_af_resistance, na.rm = TRUE),
                sum(d$priority, na.rm = TRUE), sum(d$is_novel, na.rm = TRUE)))
