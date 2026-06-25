library(tidyverse)

# Takes a SNPS data frame with snp_distribution, snp_reference, snp_position,
# and location_depth columns. Expands each distribution string into one row per
# called base, computes snp_proportion, and builds the SNP name column.
# Returns the expanded data frame (non-reference calls only, NA proportions dropped).
parse_snp_distribution <- function(snps_df) {
  snps_df$snp_distribution[is.na(snps_df$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

  snps_df <- snps_df %>% mutate(space_count = str_count(snp_distribution, " "))
  max_spaces <- max(snps_df$space_count, na.rm = TRUE)

  snps_df %>%
    relocate(snp_distribution, .after = last_col()) %>%
    separate(snp_distribution, into = paste0("Dist", 1:(1 + max_spaces)), sep = ", ", fill = "right") %>%
    pivot_longer(starts_with("Dist"), names_to = "Temp", values_to = "Dist") %>%
    select(-Temp) %>%
    filter(!is.na(Dist)) %>%
    separate(Dist, into = c("Call", "n"), sep = "=") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    mutate(SNP = paste0(snp_reference, snp_position, Call)) %>%
    filter(snp_reference != Call) %>%
    filter(!is.na(snp_proportion))
}
