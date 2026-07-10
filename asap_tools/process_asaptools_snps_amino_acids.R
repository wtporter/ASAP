#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(doParallel)
  library(foreach)
  library(parallelly)
})

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_extract_gene_table.R"))
source(file.path(.functions_dir, "_genome.snp.to.gene.snp.R"))
source(file.path(.functions_dir, "_snps.to.amino.R"))
source(file.path(.functions_dir, "_expand_codon_merges.R"))
source(file.path(.functions_dir, "_split_genbank_records.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: generate_snp_table.R <rdata> <ref1> <ref2> ...")
}

RDATA_INPUT  <- args[1]
raw_refs <- args[2:length(args)]

# RDATA_INPUT  <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/49/f254b276529b5d69d95d2f6349cb68/Combined_ASAP_Data.Rdata"
# raw_refs <- "/tgen_labs/EPIC/tporter/ASAP/nextflow/.nf-test/tests/2389292981cc5fccac0b4b3f37a2bd62/work/49/f254b276529b5d69d95d2f6349cb68/genbank_input/"

GENBANK_FILES <- c()

for (path in raw_refs) {
  if (dir.exists(path)) {
    # If the arg is a directory, get all files inside
    GENBANK_FILES <- c(GENBANK_FILES, list.files(path, full.names = TRUE, pattern = "\\.(gb|gbk|gbf|gbff|genbank)$"))
  } else if (file.exists(path)) {
    # If it's a direct file path
    GENBANK_FILES <- c(GENBANK_FILES, path)
  }
}

message("Resolved GenBank files:")
print(GENBANK_FILES)

if (length(GENBANK_FILES) == 0) {
  stop("Error: No valid GenBank files found in arguments.")
}

# Load data from ASAP_Import_XML
load(RDATA_INPUT)

SNPS <- final_snps
array_info <- final_array

# --- Initial SNPS Cleaning (as per your original code) ---
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

#Find number of SNPS
SNPS <- SNPS %>%
  mutate(space_count = str_count(snp_distribution, " "))

#Extract sub SNPs, and recalculate proportions
SNPS <- SNPS %>%
  relocate(snp_distribution, .after = last_col()) %>%
  separate(snp_distribution, into = paste0("Dist", 1:(1+max(SNPS$space_count))), sep = ", ") %>%
  pivot_longer(Dist1:ncol(.), names_to = "Temp", values_to = "Dist") %>%
  select(-Temp) %>%
  filter(!is.na(Dist)) %>%
  separate(Dist, into = c("Call", "n"), sep = "=") %>%
  mutate(snp_proportion = 100*(as.numeric(n)/as.numeric(location_depth))) %>%
  mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
  filter(snp_reference != Call) %>%
  filter(!is.na(snp_proportion))

# --- Expand codon-merged SNPs so both components get AA annotations ---
# For "complete" codon_merge pairs, replaces both individual SNPs with one
# combined-codon row (e.g. "T5118A|T5119A"); for "partial" pairs, adds the
# combined row alongside the individual rows. See _expand_codon_merges.R.
SNPS <- expand_codon_merges(SNPS, min_snp_perc = 0)

# --- Loop Through All GenBank Files ---
all_amino_acids <- list()
all_gene_snps <- list()

for (REFERENCE in GENBANK_FILES) {
  file_base <- tools::file_path_sans_ext(basename(REFERENCE))

  # Split into per-LOCUS records. genbankr::readGenBank() cannot read a
  # multi-record (multi-contig) file at all, so we never hand it the whole
  # reference; single-contig files simply yield one record.
  records <- tryCatch(split_genbank_records(REFERENCE), error = function(e) {
    message(sprintf("[WARN] Could not split %s: %s", REFERENCE, conditionMessage(e)))
    NULL
  })
  if (is.null(records)) next

  for (r in seq_len(nrow(records))) {
    locus    <- records$locus[r]
    rec_path <- records$path[r]
    # ASAP names each assay "<filebase>_<LOCUS>" (prepareJSONInput_nextflow.py),
    # so this is the exact key linking SNPs to the contig they were called on.
    assay_token <- paste0(file_base, "_", locus)

    SNPS_To_AA <- SNPS %>%
      filter(assay_name == assay_token) %>%
      select(SNP, assay_name) %>%
      distinct()

    # Only SNP-bearing contigs are parsed — critical for whole-genome refs where
    # most contigs carry no calls (e.g. one contig of a 4-contig fungal genome).
    if (nrow(SNPS_To_AA) == 0) {
      message(sprintf("No SNPs for assay %s; skipping (contig not parsed).", assay_token))
      next
    }

    message(sprintf("Processing %s (contig %s): %d SNPs", basename(REFERENCE), locus, nrow(SNPS_To_AA)))

    # Parse this contig once and share the gene table across both conversions.
    ref_df <- extract_gene_table(suppressWarnings(genbankr::readGenBank(rec_path)))

    # A contig can have SNPs but no annotated CDS (nothing to translate); skip
    # gracefully rather than crashing the downstream conversions.
    if (nrow(ref_df) == 0) {
      message(sprintf("No CDS for %s; skipping AA calling.", assay_token))
      next
    }

    # Mark genes whose reference CDS was partial (span captured pre-strip by
    # split_genbank_records). Per-gene, not per-contig: a segment may carry one
    # truncated and one complete CDS (e.g. M2 + M1).
    ref_df$partial <- paste(ref_df$start, ref_df$end) %in% records$partial_spans[[r]]
    partial_genes  <- ref_df$gene[ref_df$partial]

    gene_snps_sub <- suppressWarnings(genome.snp.to.gene.snp(
      snp_db = SNPS_To_AA, ref_seq = rec_path,
      cores = parallelly::availableCores(), ref_df = ref_df
    )) %>%
      left_join(select(SNPS_To_AA, SNP, assay_name), by = "SNP")

    amino_acids_sub <- suppressWarnings(snps.to.amino(
      snp_db = SNPS_To_AA, ref_seq = rec_path,
      cores = parallelly::availableCores(), ref_df = ref_df
    )) %>%
      left_join(select(SNPS_To_AA, SNP, assay_name), by = "SNP")

    # Flag SNPs that fall in a truncated gene so the caveat is visible in the
    # user-facing tables (Amino_Acids$AA -> "Amino Acid Change";
    # Gene_SNPS$SNP_Gene -> "SNP (Gene)"). Match by Gene, the identical string
    # in ref_df$gene and both sub-tables.
    if (length(partial_genes)) {
      note <- " (*warning: reference truncated)"
      aa_hit <- amino_acids_sub$Gene %in% partial_genes
      amino_acids_sub$AA[aa_hit] <- paste0(amino_acids_sub$AA[aa_hit], note)
      gs_hit <- gene_snps_sub$Gene %in% partial_genes
      gene_snps_sub$SNP_Gene[gs_hit] <- paste0(gene_snps_sub$SNP_Gene[gs_hit], note)
    }

    all_gene_snps[[assay_token]]   <- gene_snps_sub
    all_amino_acids[[assay_token]] <- amino_acids_sub

    message(sprintf("Success: %s", assay_token))
  }
}

# Combine results
Gene_SNPS   <- bind_rows(all_gene_snps) %>% distinct()
Amino_Acids <- bind_rows(all_amino_acids) %>% select(assay_name, SNP, Product, AA) %>% distinct()

save(Amino_Acids, Gene_SNPS, file = "SNP_Amino_Acid_Table.Rdata")
