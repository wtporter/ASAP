test_that("expand_codon_merges passes through SNPS with no codon_merge columns", {
  snps <- data.frame(
    name       = "S1",
    assay_name = "A1",
    SNP        = "C100T",
    snp_proportion = 90.0,
    stringsAsFactors = FALSE
  )
  result <- expand_codon_merges(snps, min_snp_perc = 0)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_equal(result$SNP, "C100T")
})

test_that("expand_codon_merges passes through SNPS with all-NA codon_merge_name", {
  snps <- data.frame(
    name                     = "S1",
    assay_name               = "A1",
    SNP                      = "C100T",
    snp_proportion           = 90.0,
    codon_merge_name         = NA_character_,
    codon_merge_region       = NA_character_,
    codon_merge_position     = NA_character_,
    codon_merge_codon_depth  = NA_character_,
    codon_merge_reference    = NA_character_,
    codon_merge_distribution = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- expand_codon_merges(snps, min_snp_perc = 0)

  expect_false("codon_merge_name" %in% names(result))
  expect_equal(nrow(result), 1L)
})

test_that("expand_codon_merges expands a complete codon and drops individual SNP rows", {
  snps <- make_codon_merge_snps()
  result <- expand_codon_merges(snps, min_snp_perc = 0)

  # Original row replaced by exactly one codon-derived row
  expect_equal(nrow(result), 1L)
  # codon_merge_* staging columns should be removed
  expect_false("codon_merge_name" %in% names(result))
  # TTT->TAT: only position 5118 differs, so combined name is still "T5118A"
  expect_equal(result$SNP, "T5118A")
})

test_that("expand_codon_merges drops variants below min_snp_perc", {
  snps <- make_codon_merge_snps()
  # codon_pct = 100 * 427/450 ≈ 94.9%; setting threshold above that removes it
  result <- expand_codon_merges(snps, min_snp_perc = 99)

  # No qualifying variants → original trimmed SNPS returned without codon_merge cols
  expect_false("codon_merge_name" %in% names(result))
})

test_that("make_combined_snp_name logic: single-change codon", {
  # TTT -> TAT: only position 2 changes (T->A)
  # Via expand_codon_merges on a known input
  snps <- make_codon_merge_snps()
  result <- expand_codon_merges(snps, min_snp_perc = 0)

  # TTT at positions 5117-5119, pos 2 of codon (5118) changes T->A → "T5118A"
  expect_equal(result$SNP, "T5118A")
})
