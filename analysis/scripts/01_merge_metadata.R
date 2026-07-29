#!/usr/bin/env Rscript
# 01 — assemble the genomic master table from fungiforge per-isolate outputs and
# (optionally) left-join richer sample metadata. Mirrors the bacterial 01_merge_metadata.
suppressMessages({ library(dplyr); library(readr); library(stringr) })
source(file.path(dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])), "_theme.R"))
D <- ff_dirs(); dir.create(D$output, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else file.path(D$root, "..", "results")
meta_csv    <- if (length(args) >= 2) args[[2]] else NULL   # optional external sample metadata

# gather every per-isolate master row (Stage 14 publishes *.master.tsv to 04_summary)
master_files <- Sys.glob(file.path(results_dir, "04_summary", "*.master.tsv"))
if (length(master_files) == 0)
  master_files <- Sys.glob(file.path(results_dir, "*", "14_report", "*.master.tsv"))
if (length(master_files) == 0) stop("no *.master.tsv found under ", results_dir,
                                     " — run the fungiforge pipeline first")

master <- master_files |> lapply(function(f) suppressMessages(read_tsv(f, show_col_types = FALSE))) |>
  bind_rows() |> distinct(sample, .keep_all = TRUE)
message(sprintf("[01] %d isolates from %d master files", nrow(master), length(master_files)))

if (!is.null(meta_csv) && file.exists(meta_csv)) {
  meta <- suppressMessages(read_csv(meta_csv, show_col_types = FALSE))
  key <- intersect(c("sample", "id"), names(meta))[1]
  master <- master |> left_join(meta, by = setNames(key, "sample"), suffix = c("", ".meta"))
  message(sprintf("[01] joined external metadata (%d cols) on %s", ncol(meta), key))
}

write_tsv(master, file.path(D$output, "merged_master.tsv"))
message("[01] -> output/merged_master.tsv")
