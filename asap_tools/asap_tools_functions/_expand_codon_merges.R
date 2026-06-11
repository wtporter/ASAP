expand_codon_merges <- function(SNPS) {
  library(tidyverse)

  # Nothing to do if this run produced no codon_merge annotations at all.
  if (!"codon_merge_linked_snp" %in% names(SNPS) ||
      all(is.na(SNPS$codon_merge_linked_snp) | SNPS$codon_merge_linked_snp == "")) {
    return(SNPS)
  }

  # One row per (assay_name, snp_name) node, with one sub-row per
  # <codon_merge> entry (separate_rows on the ";"-joined columns).
  merge_source <- SNPS %>%
    filter(!is.na(codon_merge_linked_snp) & codon_merge_linked_snp != "") %>%
    select(-SNP, -Call, -n, -snp_proportion, -space_count) %>%
    distinct(snp_name, assay_name, .keep_all = TRUE) %>%
    separate_rows(
      codon_merge_linked_snp, codon_merge_linkage, codon_merge_spanning_depth,
      codon_merge_variant_bases, codon_merge_variant_count, codon_merge_variant_percent,
      sep = ";"
    ) %>%
    # No <combo type="variant"> was observed for this pair -- nothing to
    # report, fall back to the individual-only annotation as if this
    # codon_merge entry didn't exist.
    filter(codon_merge_variant_bases != "")

  if (nrow(merge_source) == 0) {
    return(SNPS %>% select(-starts_with("codon_merge_")))
  }

  # SNPs where every qualifying codon_merge entry is "complete" -- their
  # individual single-position row is dropped in favor of the combined row.
  drop_keys <- merge_source %>%
    group_by(assay_name, snp_name) %>%
    summarise(all_complete = all(codon_merge_linkage == "complete"), .groups = "drop") %>%
    filter(all_complete) %>%
    transmute(assay_name, SNP = snp_name)

  # Canonical (low, high) genome-position ordering for each pair, matching
  # the <combo bases="X|Y"> convention (X = lower position, Y = higher).
  pos_linked <- as.numeric(str_extract(merge_source$codon_merge_linked_snp, "\\d+"))
  ref_linked <- str_sub(merge_source$codon_merge_linked_snp, 1, 1)
  bases      <- str_split_fixed(merge_source$codon_merge_variant_bases, "\\|", 2)
  self_is_low <- merge_source$snp_position <= pos_linked

  merge_source <- merge_source %>%
    mutate(
      mut_low      = bases[, 1],
      mut_high     = bases[, 2],
      ref_low      = if_else(self_is_low, snp_reference, ref_linked),
      pos_low      = pmin(snp_position, pos_linked),
      ref_high     = if_else(self_is_low, ref_linked, snp_reference),
      pos_high     = pmax(snp_position, pos_linked),
      combined_SNP = paste0(ref_low, pos_low, mut_low, "|", ref_high, pos_high, mut_high),
      pair_key     = paste(assay_name, pos_low, pos_high)
    )

  # Each pair appears twice (once from each SNP's perspective) -- keep one.
  pairs <- merge_source %>%
    distinct(pair_key, .keep_all = TRUE)

  extra_rows <- pairs %>%
    mutate(
      SNP             = combined_SNP,
      snp_position    = pos_low,
      snp_proportion  = as.numeric(codon_merge_variant_percent),
      snp_depth       = as.numeric(codon_merge_variant_count),
      location_depth  = as.numeric(codon_merge_spanning_depth)
    ) %>%
    select(any_of(names(SNPS)))

  SNPS %>%
    anti_join(drop_keys, by = c("assay_name", "SNP")) %>%
    bind_rows(extra_rows) %>%
    select(-starts_with("codon_merge_"))
}
