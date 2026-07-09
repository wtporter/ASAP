#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Check usage...")
}

RDATA_INPUT <- args[1]
prefix      <- args[2]
BREADTH_COVERAGE_THRESHOLD <- as.numeric(args[3])*100

# RDATA_INPUT <- "/scratch/tporter/ASAP_RSV_Results/work/1b/128f5f068cae04471a6481ab4686f7/Combined_ASAP_Data.Rdata"
# BREADTH_COVERAGE_THRESHOLD <- 80
# poi_csv     <- "NULL"
# min_depth   <- 99
# prefix      <- "ASAP_Illumina_Paired_ASAP_Tools"
# poi_csv     <- "/scratch/tporter/ASAP_TB_Validation/Updated_TB_Genes.csv"

load(RDATA_INPUT)

unique_assays <- unique(final_asap$assay_name)

any_written <- FALSE

for(ASSAY in unique_assays) {

  # Filter data for the specific assay and threshold
  Fasta_df <- final_asap %>%
    filter(assay_name == ASSAY) %>%
    filter(breadth > BREADTH_COVERAGE_THRESHOLD)

  # Check if Fasta_df is empty
  if (nrow(Fasta_df) == 0) {
    message(paste("Skipping", ASSAY, "- No samples met the breadth threshold of", BREADTH_COVERAGE_THRESHOLD))
    next # Move to the next ASSAY in the loop
  }

  # Vectorized creation of FASTA lines: much faster than a loop
  # Format: >SampleName_AssayName \n Sequence
  fasta_lines <- paste0(">", Fasta_df$name, "_", ASSAY, "\n", Fasta_df$consensus_seq)

  # Clean the assay name for a safe filename
  clean_assay <- gsub("[^[:alnum:]]", "_", ASSAY)
  file_name <- paste0("./", prefix, "_", clean_assay, ".fasta")

  # Write the file
  writeLines(fasta_lines, file_name)

  message(paste("Successfully exported", nrow(Fasta_df), "sequences to", file_name))
  any_written <- TRUE
}

if (!any_written) {
  warning_file <- paste0("./", prefix, "_WARNING_no_sequences_passed_breadth_filter.fasta")
  writeLines(
    c(
      paste0("; WARNING: No sequences were exported."),
      paste0("; No sample/assay combination exceeded the breadth coverage threshold of ",
             BREADTH_COVERAGE_THRESHOLD, "% (--breadth ", args[3], ")."),
      paste0("; Max observed breadth: ", round(max(final_asap$breadth, na.rm = TRUE), 1), "%."),
      paste0("; Lower --breadth to include more samples (e.g. --breadth 0.5 for 50%).")
    ),
    warning_file
  )
  warning(paste("No FASTA sequences exported — no samples passed breadth threshold of",
                BREADTH_COVERAGE_THRESHOLD, "%. See", warning_file))
}
