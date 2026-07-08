#!/usr/bin/env Rscript

# process_primers_to_bed.R
#
# Generate a pipeline-ready ASAP primer BED file from a 3-column primer CSV.
#
# Given a CSV with columns (primer_name, direction, sequence) and a reference
# FASTA, this locates each primer (searching both the primer and its
# reverse-complement, allowing mismatches/indels) and writes:
#   * <prefix>_primer_search_results.csv  - full search table (found + not-found,
#                                            both orientations, mismatch details)
#   * <prefix>_primers.bed                - headerless 6-column BED consumed by the
#                                            ASAP pipeline (reference, start, end,
#                                            name, score, strand); 0-based start.
#
# This script is intentionally SELF-CONTAINED (functions inlined, no source())
# so it runs identically as a standalone Rscript or from the Nextflow `bin/` PATH,
# where the shared asap_tools_functions/ helpers are not staged.
#
# Usage: process_primers_to_bed.R <primer_csv> <reference_fasta> <prefix> [max_mismatch=2]

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(foreach)
  library(data.table)
  library(doParallel)
})

###############################################################################
# Inlined helper: find.diff.in.seq
# Compare a reference string against a matched string and return the differences
# as a comma-separated "REF-POS-ALT" SNP string (or "No Difference").
###############################################################################
find.diff.in.seq <- function(a, b){

  a <- toupper(a)
  b <- toupper(b)

  seq.a <- unlist(strsplit(a, split = ""))
  seq.b <- unlist(strsplit(b, split = ""))
  diff.d <- rbind(seq.a, seq.b)
  only.diff <- diff.d[, diff.d[1, ] != diff.d[2, ]]
  pos <- which(diff.d[1, ] != diff.d[2, ])

  Out <- as.data.frame(t(rbind(pos, as.data.frame(only.diff))))

  if(nrow(Out) == 0)
  {Out = "No Difference"} else {
    Out <- mutate(Out, SNP = paste(seq.a, `1`, seq.b, sep = "-"))
    Out <- paste(Out$SNP, collapse = ", ")
  }

  return(Out)
}

###############################################################################
# Inlined helper: find.primers
# Search a reference FASTA for each primer (forward + reverse-complement),
# returning a data frame with match coordinates and per-match mismatch details.
###############################################################################
find.primers <- function(fasta.path, primer.names, primer.direction, primer.list, max.mismatch, cores = 1){

  # Read the FASTA file
  sequence_file <- Biostrings::readDNAStringSet(fasta.path)

  # Search primers in parallel when more than one core is available. On Unix,
  # registerDoParallel(cores) uses forked workers, so the reference is shared
  # copy-on-write (no per-worker reload). Fall back to sequential otherwise.
  cores <- max(1L, min(as.integer(cores), length(primer.list)))
  if (cores > 1L) {
    doParallel::registerDoParallel(cores = cores)
    `%run%` <- foreach::`%dopar%`
  } else {
    `%run%` <- foreach::`%do%`
  }
  on.exit(if (cores > 1L) doParallel::stopImplicitCluster(), add = TRUE)

  # Iterate over each primer
  bed_data <- foreach(i = 1:length(primer.list), .combine = rbind,
                      .packages = "Biostrings") %run% {

    PRIMER <- primer.names[[i]]
    PRIMER_DIRECTION <- primer.direction[[i]]
    PRIMER_SEQUENCE <- primer.list[[i]]

    Combine_Out <- data.frame()

    for (GENE in 1:length(sequence_file)) {

      sequences <- sequence_file[[GENE]]

      ASSAY <- names(sequence_file)[[GENE]]

      # fixed = "subject": treat IUPAC ambiguity codes in the PRIMER as wildcards
      # while keeping the reference literal. This requires a DNAString subject
      # (passing as.character() disables ambiguity handling and errors on
      # fixed="subject"), so the raw DNAString `sequences` is passed directly.
      Temp <- matchPattern(pattern = PRIMER_SEQUENCE, subject = sequences,
                           max.mismatch = max.mismatch, with.indels = T, fixed = "subject",
                           algorithm = "auto")

      if(nrow(data.frame(Temp)) > 0) {

        Out <- cbind(assay = ASSAY,
                     primer = PRIMER,
                     direction = PRIMER_DIRECTION,
                     primer.sequence = PRIMER_SEQUENCE,
                     primer.search = PRIMER_SEQUENCE,
                     searched.direction = "Forward",
                     data.frame(Temp))

      }else{
        Out <- cbind(assay = ASSAY,
                     primer = PRIMER,
                     direction = PRIMER_DIRECTION,
                     primer.sequence = PRIMER_SEQUENCE,
                     primer.search = PRIMER_SEQUENCE,
                     searched.direction = "Forward",
                     start = NA,
                     end = NA,
                     width = NA,
                     seq = NA)
      }

      ## Do reverse compliment search
      Temp <- matchPattern(pattern = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))), subject = sequences,
                           max.mismatch = max.mismatch, with.indels = T, fixed = "subject",
                           algorithm = "auto")

      if(nrow(data.frame(Temp)) > 0) {

        Out <- rbind(Out, cbind(assay = ASSAY,
                                primer = PRIMER,
                                direction = PRIMER_DIRECTION,
                                primer.sequence = PRIMER_SEQUENCE,
                                primer.search = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))),
                                searched.direction = "Reverse",
                                data.frame(Temp)))

      }else{
        Out <- rbind(Out, cbind(assay = ASSAY,
                                primer = PRIMER,
                                direction = PRIMER_DIRECTION,
                                primer.sequence = PRIMER_SEQUENCE,
                                primer.search = as.character(reverseComplement(DNAString(PRIMER_SEQUENCE))),
                                searched.direction = "Reverse",
                                start = NA,
                                end = NA,
                                width = NA,
                                seq = NA))
      }
      Combine_Out <- rbind(Combine_Out, Out)
    }

    Combine_Out
  }

  bed_data$Mismatches <- NA

  for (ROW in which(complete.cases(bed_data$seq))) {
    bed_data$Mismatches[ROW] <- find.diff.in.seq(bed_data$primer.search[ROW], bed_data$seq[ROW])
  }

  return(bed_data)
}

###############################################################################
# Inlined + adapted helper: create.asap.bed.file
# Emit a headerless 6-column standard BED (reference, start, end, name, score,
# strand) as required by the ASAP pipeline (maskPrimers.py / iVar / SNP table).
# Coordinates are converted from Biostrings' 1-based match start to 0-based BED
# start, and strand is derived from the primer direction (F -> +, R -> -).
###############################################################################
create.asap.bed.file <- function(find.primers.df, save.path = NA){

  # Keep only found primers (drops the NA rows for orientations/records with no match)
  Temp <- find.primers.df %>%
    na.omit()

  # Build the 6 standard BED columns
  Temp <- Temp %>%
    transmute(
      reference = assay,
      start     = as.integer(start) - 1L,   # 1-based match start -> 0-based BED start
      end       = as.integer(end),
      name      = primer,
      score     = 0L,
      strand    = dplyr::recode(as.character(direction),
                                "F" = "+", "f" = "+",
                                "R" = "-", "r" = "-",
                                .default = as.character(direction))
    )

  # Write a headerless, tab-delimited BED
  if(!is.na(save.path)){
    write.table(Temp, file = save.path, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  return(Temp)
}

###############################################################################
# --- Argument parsing ---
###############################################################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: process_primers_to_bed.R <primer_csv> <reference_fasta> <prefix> [max_mismatch=2] [cores=1]")
}

PRIMER_CSV <- args[1]
REF_FASTA  <- args[2]
PREFIX     <- args[3]
MAX_MM     <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 2
CORES      <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1L

###############################################################################
# --- Read + validate the primer CSV ---
###############################################################################
primers <- readr::read_csv(PRIMER_CSV, show_col_types = FALSE)

required_cols <- c("primer_name", "direction", "sequence")
missing_cols  <- setdiff(required_cols, colnames(primers))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Primer CSV '%s' is missing required column(s): %s. Expected columns: primer_name, direction, sequence.",
    PRIMER_CSV, paste(missing_cols, collapse = ", ")))
}

###############################################################################
# --- Run the primer search and build the BED ---
###############################################################################
search_results <- find.primers(
  fasta.path       = REF_FASTA,
  primer.names     = primers$primer_name,
  primer.direction = primers$direction,
  primer.list      = primers$sequence,
  max.mismatch     = MAX_MM,
  cores            = CORES
)

results_path <- paste0(PREFIX, "_primer_search_results.csv")
bed_path     <- paste0(PREFIX, "_primers.bed")
summary_path <- paste0(PREFIX, "_primer_match_summary.csv")

# Full search table (found + not-found, both orientations) for QC
readr::write_csv(search_results, results_path)

# Pipeline-ready 6-column headerless BED
bed <- create.asap.bed.file(search_results, save.path = bed_path)

if (nrow(bed) == 0) {
  warning(sprintf(
    "No primers were located in '%s' (max.mismatch=%s). The BED '%s' is empty.",
    REF_FASTA, MAX_MM, bed_path))
}

###############################################################################
# --- Match summary report: matches per primer per reference (wide) ---
# One row per primer (in input order), one count column per reference sequence,
# plus a total_matches column. Primers that matched nothing show an all-zero row.
###############################################################################
ref_names <- names(Biostrings::readDNAStringSet(REF_FASTA))

# Count found matches per (primer, direction, reference)
match_counts <- search_results %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::count(primer_name = primer, direction = direction, reference = assay,
               name = "n")

# Base = every input primer, preserving input order
base <- primers %>%
  dplyr::transmute(primer_name, direction) %>%
  dplyr::mutate(.order = dplyr::row_number())

match_summary <- match_counts %>%
  tidyr::pivot_wider(names_from = reference, values_from = n, values_fill = 0)

match_summary <- base %>%
  dplyr::left_join(match_summary, by = c("primer_name", "direction"))

# Ensure a column exists for every reference in the FASTA (even 0-match refs),
# fill missing counts with 0, and order rows by original input order.
for (rn in ref_names) if (!rn %in% colnames(match_summary)) match_summary[[rn]] <- 0L
match_summary <- match_summary %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(ref_names), ~tidyr::replace_na(as.integer(.), 0L))) %>%
  dplyr::arrange(.order) %>%
  dplyr::mutate(total_matches = rowSums(dplyr::across(dplyr::all_of(ref_names)))) %>%
  dplyr::select(primer_name, direction, dplyr::all_of(ref_names), total_matches)

readr::write_csv(match_summary, summary_path)

cat(sprintf("Wrote %d primer BED entries to %s\n", nrow(bed), bed_path))
cat(sprintf("Wrote primer search results to %s\n", results_path))
cat(sprintf("Wrote primer match summary to %s\n", summary_path))
