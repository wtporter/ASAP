#!/usr/bin/env Rscript

# combine_masked_reads_wide.R
#
# Combine the per-sample masked-reads-per-primer TSVs into a single WIDE table:
# one row per (ref_name, primer_name, direction) and one column per sample, where
# each cell is the number of reads masked for that primer in that sample.
#
# Each per-sample input has a long-format header:
#   sample_id  ref_name  primer_name  direction  masked_reads
# so the sample identity travels with the data and no filename parsing is needed.
#
# Sample columns are ordered alphabetically for deterministic, diff-friendly output.
# Missing (primer x sample) combinations are filled with 0.
#
# Usage: combine_masked_reads_wide.R <output_tsv> <per_sample_tsv> [<per_sample_tsv> ...]

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: combine_masked_reads_wide.R <output_tsv> <per_sample_tsv> [<per_sample_tsv> ...]")
}

out_path <- args[[1]]
in_files <- args[-1]

id_cols <- c("ref_name", "primer_name", "direction")

# Read + stack every per-sample long-format table. Read all columns as character
# first so an empty/edge-case file can't force a bad column type, then coerce.
long <- purrr::map_dfr(in_files, ~ readr::read_tsv(
  .x, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character())
))

# Guard: nothing to combine (e.g. no primers matched any sample) -> write a
# header-only file with just the id columns so downstream steps still find it.
if (nrow(long) == 0) {
  readr::write_tsv(tibble::tibble(!!!setNames(rep(list(character()), length(id_cols)), id_cols)),
                   out_path)
  cat(sprintf("No masked-reads records found; wrote empty wide table to %s\n", out_path))
  quit(save = "no", status = 0)
}

long <- long %>%
  dplyr::mutate(masked_reads = as.integer(masked_reads))

wide <- long %>%
  tidyr::pivot_wider(
    id_cols     = dplyr::all_of(id_cols),
    names_from  = sample_id,
    values_from = masked_reads,
    values_fill = 0L
  )

# Deterministic column + row ordering
sample_cols <- sort(setdiff(colnames(wide), id_cols))
wide <- wide %>%
  dplyr::select(dplyr::all_of(id_cols), dplyr::all_of(sample_cols)) %>%
  dplyr::arrange(ref_name, primer_name, direction)

readr::write_tsv(wide, out_path)

cat(sprintf("Wrote wide masked-reads table (%d primers x %d samples) to %s\n",
            nrow(wide), length(sample_cols), out_path))
