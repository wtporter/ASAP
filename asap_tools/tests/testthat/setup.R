library(testthat)
library(tidyverse)
library(Biostrings)

# Functions dir: override with ASAP_FUNCTIONS_DIR env var, otherwise derive from cwd
# cwd when run via test_dir() is tests/testthat/; go up two levels to asap_tools/
.functions_dir <- Sys.getenv("ASAP_FUNCTIONS_DIR",
  unset = file.path(dirname(dirname(getwd())), "asap_tools_functions"))

for (f in c("_snps.to.amino.R", "_genome.snp.to.gene.snp.R", "_expand_codon_merges.R",
            "_read.ASAP.individual.R", "_read.ASAP.snps.individual.R",
            "_parse_snp_distribution.R")) {
  source(file.path(.functions_dir, f))
}

.fixture_dir <- file.path(getwd(), "fixtures")

# ---------------------------------------------------------------------------
# Inline fixture builders — no external files needed for pure unit tests
# ---------------------------------------------------------------------------

make_snps_df <- function() {
  data.frame(
    name             = "TEST_SAMPLE",
    assay_name       = "TB_rpoB",
    snp_position     = c(1349L, 1547L),
    snp_reference    = c("C", "A"),
    location_depth   = c(450L, 448L),
    snp_distribution = c("A=0, C=45, G=0, T=405, _=0", "A=40, C=0, G=0, T=408, _=0"),
    stringsAsFactors = FALSE
  )
}

make_codon_merge_snps <- function() {
  data.frame(
    name                     = "TEST_SAMPLE",
    assay_name               = "TB_rpoB",
    SNP                      = "T5118A",
    snp_position             = 5118L,
    snp_proportion           = 95.0,
    snp_depth                = 427L,
    location_depth           = 450L,
    codon_merge_name         = "rpoB_codon_1706",
    codon_merge_region       = "5117-5119",
    codon_merge_position     = "5117-5119",
    codon_merge_codon_depth  = "450",
    codon_merge_reference    = "TTT",
    codon_merge_distribution = "TTT=23,TAT=427",
    stringsAsFactors = FALSE
  )
}
