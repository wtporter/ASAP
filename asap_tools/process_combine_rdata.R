#!/usr/bin/env Rscript

# Load necessary libraries
# tidyverse for data manipulation, data.table for high-speed binding
suppressPackageStartupMessages({
  library(tidyverse)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(parallelly)
})

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_shorten_sample_names.R"))

# 1. Capture Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: process_combine_rdata.R <optional_poi_csv> <file_name> <input_files...>", call. = FALSE)
}

poi_csv   <- args[1]
file_name <- args[2]
files     <- args[3:length(args)]


# 2. Setup Parallel Backend
# parallelly::availableCores() is SLURM-aware and respects cpus allocated to the job
num_cores <- parallelly::availableCores()
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# Only use parallel if more than 1 file
if(length(files) > 1) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
} else {
  registerDoSEQ() # Fallback to sequential for single file
}

message(paste("🚀 Combining", length(files), "samples using", num_cores, "cores..."))

# 3. Parallel Loading with foreach
# We use %dopar% to read files simultaneously across available cores
combined_list <- foreach(f = files, .packages = c("tidyverse")) %dopar% {
  
  if (!file.exists(f)) {
    stop(paste("File not found in task directory:", f))
  }
  
  # Create a local environment for each file to prevent object collision
  temp_env <- new.env()
  load(f, envir = temp_env)
  
  # 2. Handle Positions of Interest (Optional)
  if (is.na(poi_csv) || poi_csv == "NULL" || poi_csv == "" || is.null(poi_csv)) {
    message("No Positions of Interest provided. Generating reference data...")

    } else {
    message(paste("Loading positions of interest from:", poi_csv))
    
    if (!file.exists(poi_csv)) {
      stop(paste("Positions of interest not found:", poi_csv))
    }
      
    genes <- read.csv(poi_csv)
    
    if (nrow(genes) == 0) {
      stop(paste("Positions of interest file is seemingly empty. Please check:", poi_csv))
    }
    
    # Generate positions for each gene in the CSV
    Gene_Positions <- genes %>%
      rowwise() %>%
      do(data.frame(
        position = seq(min(.$start, .$end), max(.$start, .$end)), 
        gene = .$gene,
        reference = .$seqnames
      )) %>%
      ungroup()
    
    # Guard: the POI CSV's `seqnames` must match the data's `assay_name`.
    # If they don't overlap at all, the semi_join below silently drops every
    # row, leaving final_array empty and crashing downstream coverage steps
    # ~50 minutes later. Fail loudly here, naming both sides of the mismatch.
    poi_refs   <- unique(as.character(Gene_Positions$reference))
    data_names <- unique(as.character(temp_env$array_info$assay_name))
    if (length(data_names) > 0 && !any(data_names %in% poi_refs)) {
      stop(paste0(
        "Positions-of-interest reference names do not match the data.\n",
        "  POI CSV 'seqnames': ", paste(poi_refs,   collapse = ", "), "\n",
        "  data 'assay_name':  ", paste(data_names, collapse = ", "), "\n",
        "Fix the 'seqnames' column in ", poi_csv,
        " to match the reference/assay name."
      ))
    }

    # filter unneeded array info...
    temp_env$array_info <- temp_env$array_info %>%
      semi_join(Gene_Positions, by = c("position" = "position", "assay_name" = "reference"))
    
    # Remove array_info fields from ASAP
    temp_env$ASAP$n_reads <- NULL
    temp_env$ASAP$quality_discards <- NULL
    temp_env$ASAP$proportions <- NULL
    temp_env$ASAP$depths <- NULL
    
    # Remove unneeded data from fields.
    # For TB Example, no way to easily clean SNPs to match POI filtering... 
  }
  
  # Return as a structured list for easier extraction
  list(
    asap = temp_env$ASAP,
    snps = temp_env$SNPS,
    info = temp_env$array_info
  )
}

# Explicitly stop the cluster to free system resources
if (length(files) > 1) stopCluster(cl)

# 4. Fast Binding with data.table
# rbindlist is written in C and is significantly faster than map_df or rbind
message("📊 Merging data frames...")

# Merge ASAP summary data
final_asap  <- as.data.frame(data.table::rbindlist(map(combined_list, "asap"), fill = TRUE))
gc() # Clear memory immediately after large merges

# Merge SNP data
final_snps  <- as.data.frame(data.table::rbindlist(map(combined_list, "snps"), fill = TRUE))
gc()

# Merge large Array Info (Depth/Proportions)
final_array <- as.data.frame(data.table::rbindlist(map(combined_list, "info"), fill = TRUE))
gc()

# Compute shortened plot labels once against the full cohort so every
# downstream script/plot uses the same name -> name_short mapping. The full
# `name` column is left untouched (still used for exports/tables/traceability).
name_map <- shorten_sample_names(final_asap$name)
final_asap$name_short  <- name_map[final_asap$name]
final_snps$name_short  <- name_map[final_snps$name]
final_array$name_short <- name_map[final_array$name]

# 5. Save Combined Outputs
# Saving both as a compressed Rdata object and a flat CSV summary
message(paste("💾 Saving results"))

save(final_asap, final_snps, final_array, file = paste0(file_name, "_ASAP_Data.Rdata"))
cols_to_drop <- c("consensus_seq", "depths", "proportions", "quality_discards", "n_reads")
write.csv(final_asap[, !names(final_asap) %in% cols_to_drop, drop = FALSE], file = paste0(file_name, "_Summary.csv"), row.names = FALSE)

message("✅ Success: Combined data saved to current working directory.")