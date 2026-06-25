library(tidyverse)

expand_codon_merges <- function(SNPS, min_snp_perc = 0) {
  # Nothing to do if this run produced no codon_merge annotations at all.
  if (!"codon_merge_name" %in% names(SNPS) ||
      all(is.na(SNPS$codon_merge_name) | SNPS$codon_merge_name == "")) {
    return(SNPS %>% select(-any_of(starts_with_codon_merge_cols(SNPS))))
  }

  # ---------------------------------------------------------------------------
  # A. Expand semicolons: one row per (SNP node, codon_merge entry)
  # ---------------------------------------------------------------------------
  merge_source <- SNPS %>%
    filter(!is.na(codon_merge_name) & codon_merge_name != "") %>%
    separate_rows(
      codon_merge_name, codon_merge_region, codon_merge_position,
      codon_merge_codon_depth, codon_merge_reference, codon_merge_distribution,
      sep = ";"
    ) %>%
    filter(!is.na(codon_merge_distribution) & codon_merge_distribution != "")

  if (nrow(merge_source) == 0) {
    return(SNPS %>% select(-starts_with("codon_merge_")))
  }

  # ---------------------------------------------------------------------------
  # B. Deduplicate by codon: the same codon_merge appears on every SNP in
  #    the codon (e.g. T5118A and T5119A both carry ORF1ab_codon_1618).
  #    Keep one representative row per (assay_name, sample name, codon).
  # ---------------------------------------------------------------------------
  codon_rows <- merge_source %>%
    distinct(assay_name, name, codon_merge_name, .keep_all = TRUE)

  # ---------------------------------------------------------------------------
  # C & D. Parse distribution, filter non-reference codons above threshold
  # ---------------------------------------------------------------------------
  codon_variants <- codon_rows %>%
    separate_rows(codon_merge_distribution, sep = ",") %>%
    separate(codon_merge_distribution,
             into = c("codon_seq", "codon_count"),
             sep  = "=", extra = "merge", fill = "right") %>%
    mutate(
      codon_count     = as.numeric(codon_count),
      codon_depth_num = as.numeric(codon_merge_codon_depth),
      codon_pct       = 100 * codon_count / pmax(codon_depth_num, 1)
    ) %>%
    filter(codon_seq != codon_merge_reference) %>%
    filter(codon_pct >= min_snp_perc)

  if (nrow(codon_variants) == 0) {
    return(SNPS %>% select(-starts_with("codon_merge_")))
  }

  # ---------------------------------------------------------------------------
  # E. Build combined SNP name from position range + reference + variant codon
  #    e.g. position="5117-5119", reference="TTT", variant="TAA" → "T5118A|T5119A"
  # ---------------------------------------------------------------------------
  make_combined_snp_name <- function(position_str, ref_codon, var_codon) {
    bounds    <- as.numeric(strsplit(position_str, "-")[[1]])
    positions <- seq(bounds[1], bounds[2])
    ref_bases <- strsplit(ref_codon, "")[[1]]
    var_bases <- strsplit(var_codon, "")[[1]]
    changed   <- which(ref_bases != var_bases)
    if (length(changed) == 0) return(paste0(ref_codon, bounds[1], var_codon))
    paste(paste0(ref_bases[changed], positions[changed], var_bases[changed]),
          collapse = "|")
  }

  codon_variants <- codon_variants %>%
    rowwise() %>%
    mutate(combined_SNP = make_combined_snp_name(
      codon_merge_position, codon_merge_reference, codon_seq
    )) %>%
    ungroup()

  # ---------------------------------------------------------------------------
  # F. Determine linkage per (assay_name, sample, codon):
  #    complete = exactly one non-reference codon is above threshold
  # ---------------------------------------------------------------------------
  linkage <- codon_variants %>%
    group_by(assay_name, name, codon_merge_name) %>%
    summarise(n_variants = n(), .groups = "drop") %>%
    mutate(is_complete = n_variants == 1)

  codon_variants <- codon_variants %>%
    left_join(linkage %>% select(assay_name, name, codon_merge_name, is_complete),
              by = c("assay_name", "name", "codon_merge_name"))

  # ---------------------------------------------------------------------------
  # G. For complete codons: identify and drop the individual per-base SNP rows
  #    that belong to those codons (every SNP node that references this codon).
  # ---------------------------------------------------------------------------
  drop_keys <- merge_source %>%
    inner_join(linkage %>% filter(is_complete) %>%
                 select(assay_name, name, codon_merge_name),
               by = c("assay_name", "name", "codon_merge_name")) %>%
    distinct(assay_name, name, SNP)

  SNPS_trimmed <- SNPS %>%
    anti_join(drop_keys, by = c("assay_name", "name", "SNP"))

  # ---------------------------------------------------------------------------
  # H. Build extra rows — one per qualifying non-reference codon variant
  # ---------------------------------------------------------------------------
  extra_rows <- codon_variants %>%
    mutate(
      SNP            = combined_SNP,
      snp_position   = as.numeric(str_extract(codon_merge_position, "^\\d+")),
      snp_proportion = codon_pct,
      snp_depth      = codon_count,
      location_depth = codon_depth_num
    ) %>%
    select(any_of(names(SNPS)))

  # ---------------------------------------------------------------------------
  # I. Bind and drop all codon_merge_* staging columns
  # ---------------------------------------------------------------------------
  bind_rows(SNPS_trimmed, extra_rows) %>%
    select(-starts_with("codon_merge_"))
}

# Helper: column names starting with "codon_merge_" that exist in a data frame
starts_with_codon_merge_cols <- function(df) {
  names(df)[startsWith(names(df), "codon_merge_")]
}
