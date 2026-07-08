library(tidyverse)

# Takes a SNPS data frame with snp_distribution, snp_reference, snp_position,
# and location_depth columns. Expands each distribution string into one row per
# called base, computes snp_proportion, and builds the SNP name column.
# Returns the expanded data frame (non-reference calls only, NA proportions dropped).
parse_snp_distribution <- function(snps_df) {
  snps_df$snp_distribution[is.na(snps_df$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

  # separate_longer_delim expands each distribution to one row per allele directly.
  # Avoid separate()+pivot_longer(): that first builds a (1 + max spaces)-wide frame
  # for every row, exploding to a transient ~150M-row / >20GB peak when indel-rich
  # positions push the allele count to 50+.
  snps_df %>%
    separate_longer_delim(snp_distribution, delim = ", ") %>%
    separate_wider_delim(snp_distribution, delim = "=", names = c("Call", "n"),
                         too_many = "merge", too_few = "align_start") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
    filter(snp_reference != Call) %>%
    filter(!is.na(snp_proportion))
}
