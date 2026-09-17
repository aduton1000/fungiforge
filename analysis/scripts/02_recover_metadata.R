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

# v0.2 columns that an older results directory will not have
for (col in c("mycotoxin_clusters", "bioactive_clusters", "te_percent", "n_mycovirus", "n_secreted", "n_effectors",
              "n_cazymes", "n_phibase_hits", "mito_size_kb", "resistance_read_support", "gate", "species_confidence"))
  if (!col %in% names(m)) m[[col]] <- NA

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
         priority = str_detect(species, regex("Candida auris|Aspergillus fumigatus|Cryptococcus|Candida albicans|Candida glabrata|Nakaseomyces|Fusarium", ignore_case = TRUE)),
         # W3.3: columns the pipeline populates from v0.2 onwards (absent in older tables -> NA)
         gate_passed = !str_detect(coalesce(gate, "pass"), "^skipped"),
         id_high_conf = coalesce(species_confidence, "") == "high",
         mycotoxin_clusters = ifelse(is.na(mycotoxin_clusters), "NA", mycotoxin_clusters),
         has_mycotoxin = !mycotoxin_clusters %in% c("none", "NA"),
         n_mycotoxin_clusters = ifelse(has_mycotoxin, str_count(mycotoxin_clusters, ";") + 1L, 0L),
         aflatoxin = str_detect(mycotoxin_clusters, "aflatoxin"),
         gliotoxin = str_detect(mycotoxin_clusters, "gliotoxin"),
         te_percent = suppressWarnings(as.numeric(te_percent)),
         n_mycovirus = suppressWarnings(as.integer(n_mycovirus)),
         n_secreted = suppressWarnings(as.integer(n_secreted)),
         n_effectors = suppressWarnings(as.integer(n_effectors)),
         n_cazymes = suppressWarnings(as.integer(n_cazymes)),
         n_phibase_hits = suppressWarnings(as.integer(n_phibase_hits)),
         mito_size_kb = suppressWarnings(as.numeric(mito_size_kb)),
         resistance_read_confirmed = str_detect(coalesce(resistance_read_support, ""), "confirmed:[1-9]"))

write_csv(d, file.path(D$output, "analysis_dataset.csv"))
message(sprintf("[02] analysis_dataset.csv: %d isolates | %d compartment-known | %d AF-resistant | %d priority | %d novel | %d mycotoxigenic | %d gated out",
                nrow(d), sum(!is.na(d$compartment)), sum(d$has_af_resistance, na.rm = TRUE),
                sum(d$priority, na.rm = TRUE), sum(d$is_novel, na.rm = TRUE),
                sum(d$has_mycotoxin, na.rm = TRUE), sum(!d$gate_passed, na.rm = TRUE)))
