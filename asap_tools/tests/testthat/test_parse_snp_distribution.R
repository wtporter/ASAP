test_that("parse_snp_distribution expands distribution strings into rows", {
  snps <- make_snps_df()
  result <- parse_snp_distribution(snps)

  expect_s3_class(result, "data.frame")
  # Each row had 4 non-reference alleles possible (A,C,G,T,_ minus ref);
  # we filter out ref==Call and NA proportions, so should have multiple rows
  expect_gt(nrow(result), nrow(snps))
})

test_that("parse_snp_distribution adds snp_proportion column", {
  result <- parse_snp_distribution(make_snps_df())

  expect_true("snp_proportion" %in% names(result))
  expect_true(is.numeric(result$snp_proportion))
  expect_true(all(result$snp_proportion >= 0 & result$snp_proportion <= 100))
})

test_that("parse_snp_distribution adds SNP name column in <ref><pos><call> format", {
  result <- parse_snp_distribution(make_snps_df())

  expect_true("SNP" %in% names(result))
  expect_true(all(grepl("^[ACGT_]\\d+[ACGT_]$", result$SNP)))
})

test_that("parse_snp_distribution fills NA snp_distribution with zeros", {
  snps <- make_snps_df()
  snps$snp_distribution[1] <- NA
  result <- parse_snp_distribution(snps)

  # Row 1 had NA distribution, replaced by "A=0, T=0, C=0, G=0, _=0"
  # All resulting proportions for that row should be 0
  row1_results <- result[result$snp_position == snps$snp_position[1], ]
  expect_true(all(row1_results$snp_proportion == 0))
})

test_that("parse_snp_distribution computes correct proportion for known input", {
  snps <- make_snps_df()
  # Row 1: C1349, depth 450, "A=0, C=45, G=0, T=405, _=0"
  # T call: 405/450 * 100 = 90.0
  result <- parse_snp_distribution(snps)

  t_row <- result[result$snp_position == 1349L & result$Call == "T", ]
  expect_equal(nrow(t_row), 1L)
  expect_equal(round(t_row$snp_proportion, 1), 90.0)
})

test_that("parse_snp_distribution filters out reference==call rows", {
  snps <- make_snps_df()
  # snp_reference[1] = "C": rows where Call == "C" should be dropped
  result <- parse_snp_distribution(snps)

  ref_rows <- result[result$snp_position == 1349L & result$Call == "C", ]
  expect_equal(nrow(ref_rows), 0L)
})
