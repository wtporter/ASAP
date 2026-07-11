#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: generate_coverage_table.R <combined_rdata> <min_depth> <prefix> <optional_poi_csv>")
}

rdata_input <- args[1]
min_depth   <- as.numeric(args[2])
prefix      <- args[3]
poi_csv     <- args[4] # Will be "NULL" if not provided

poi_csv <- if(length(args) >= 4) args[4] else "NULL"

# rdata_input <- "/scratch/tporter/ASAP_SC2_Validation/ASAP_Illumina_Paired_SE_Test/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "NULL"
# rdata_input <- "/scratch/tporter/ASAP_TB_Validation/ASAP_TB_Subset/ASAP_R_Data/Combined_ASAP_Data.Rdata"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "/scratch/tporter/ASAP_TB_Validation/Updated_TB_Genes.csv"


# 1. Load Data
load(rdata_input) # Loads final_asap, final_snps, final_array

# 2. Handle Positions of Interest (Optional)
if (is.na(poi_csv) || poi_csv == "NULL" || poi_csv == "") {
  message("No Positions of Interest provided. Generating Whole-Reference coverage summary.")

  # Denominator = exact reference length from the full-length consensus sequence in final_asap.
  # The consensus stays full-length even when --prune-per-base sparsifies the numeric arrays, so
  # counting final_array rows (n()) would otherwise collapse Coverage to ~100% for pruned runs.
  ref_lengths <- final_asap %>%
    filter(consensus_seq != "No Consensus Sequence") %>%
    group_by(name, assay_name) %>%
    summarise(total_bp = sum(nchar(consensus_seq)), .groups = 'drop')

  Amplicon_Coverage <- final_array %>%
    group_by(name, assay_name) %>%
    summarise(n_cov = sum(depth >= min_depth, na.rm = TRUE), .groups = 'drop') %>%
    left_join(ref_lengths, by = c("name", "assay_name")) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, Coverage) %>%
    pivot_wider(names_from = assay_name, values_from = Coverage)


} else {
  message(paste("Loading positions of interest from:", poi_csv))
  genes <- read.csv(poi_csv)

  # Generate positions for each gene in the CSV
  Gene_Positions <- genes %>%
    rowwise() %>%
    do(data.frame(
      position = seq(min(.$start, .$end), max(.$start, .$end)),
      gene = .$gene,
      assay_name = .$seqnames
    )) %>%
    ungroup()

  # Guard: POI `seqnames` must match the data's `assay_name`, otherwise the
  # left_join + filter below drops every row and the styling loop crashes on an
  # empty table. Fail loudly with both name sets (mirrors process_combine_rdata.R).
  poi_refs   <- unique(as.character(Gene_Positions$assay_name))
  data_names <- unique(as.character(final_array$assay_name))
  if (length(data_names) > 0 && !any(data_names %in% poi_refs)) {
    stop(paste0(
      "Positions-of-interest reference names do not match the data.\n",
      "  POI CSV 'seqnames': ", paste(poi_refs,   collapse = ", "), "\n",
      "  data 'assay_name':  ", paste(data_names, collapse = ", "), "\n",
      "Fix the 'seqnames' column in ", poi_csv,
      " to match the reference/assay name."
    ))
  }

  # Build a full (sample x gene-position) grid so every gene position is counted in the denominator
  # even when --prune-per-base drops uncovered positions from final_array (n() over the sparse array
  # would otherwise undercount total_bp and inflate Coverage). Scoped to the sample/assay pairs
  # actually present so samples that lack an assay aren't invented; depth is NA where the sample has
  # no coverage at a gene position, so it correctly counts as uncovered.
  sample_assays <- final_array %>% distinct(run, name, assay_name)
  array_info <- sample_assays %>%
    left_join(Gene_Positions, by = "assay_name", relationship = "many-to-many") %>%
    filter(!is.na(gene)) %>%
    left_join(final_array, by = c("run", "name", "assay_name", "position"))

  Assay_Coverage <- array_info %>%
    group_by(name, assay_name) %>%
    summarise(
      total_bp = n(),
      n_cov = sum(depth >= min_depth, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, Coverage) %>%
    pivot_wider(names_from = assay_name, values_from = Coverage)

  POI_Coverage <- array_info %>%
    group_by(name, assay_name, gene) %>%
    summarise(
      total_bp = n(),
      n_cov = sum(depth >= min_depth, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Coverage = round(100 * (n_cov / total_bp), 2)) %>%
    select(name, assay_name, gene, Coverage) %>%
    pivot_wider(names_from = c(assay_name, gene), values_from = Coverage)

  Amplicon_Coverage <- full_join(Assay_Coverage, POI_Coverage)

}

# 5. Styling and Excel Export
wb <- createWorkbook("TGen North")
addWorksheet(wb, "Amplicon_Coverage")

color_breaks_10 <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850", "#006837")

getStyle_simple_100 <- function(value) {
  # Define the border vector once
  all_borders <- c("top", "bottom", "left", "right")
  if (is.na(value)) {
    return(createStyle(fgFill = "gray75", border = all_borders))
  }
  clamped_value <- max(0, min(100, value))
  color_index <- max(1, min(10, ceiling((clamped_value + 1e-6) / 10)))
  return(createStyle(fgFill = color_breaks_10[color_index], border = all_borders))
}

# Apply styles (Starting from col 2 to skip sample name). Guard against an empty
# or single-column table so `1:nrow`/`2:ncol` can't run past the end — use
# seq_len/seq so a zero-row or no-coverage-column result is skipped, not crashed.
if (nrow(Amplicon_Coverage) == 0 || ncol(Amplicon_Coverage) < 2) {
  warning("Amplicon_Coverage has no coverage values to style ",
          "(rows: ", nrow(Amplicon_Coverage), ", cols: ", ncol(Amplicon_Coverage),
          "). Writing the table without conditional styling.")
} else {
  for (row in seq_len(nrow(Amplicon_Coverage))) {
    for (col in seq(2, ncol(Amplicon_Coverage))) {
      val <- Amplicon_Coverage[[row, col]]
      addStyle(wb, "Amplicon_Coverage", style = getStyle_simple_100(val), rows = row + 1, cols = col)
    }
  }
}

writeData(wb, "Amplicon_Coverage", Amplicon_Coverage)
saveWorkbook(wb, paste0(prefix, "_Coverage_Report.xlsx"), overwrite = TRUE)
