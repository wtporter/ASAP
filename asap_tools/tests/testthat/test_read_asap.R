test_that("read.ASAP.individual returns expected columns from fixture XML", {
  xml_path <- file.path(.fixture_dir, "sample_asap.xml")
  skip_if_not(file.exists(xml_path), "fixture XML not found")

  result <- suppressMessages(read.ASAP.individual(xml_path))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0L)

  expected_cols <- c("name", "assay_name", "breadth", "avg_depth",
                     "depths", "proportions", "total_reads")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("read.ASAP.individual extracts correct sample name", {
  xml_path <- file.path(.fixture_dir, "sample_asap.xml")
  skip_if_not(file.exists(xml_path), "fixture XML not found")

  result <- suppressMessages(read.ASAP.individual(xml_path))

  expect_equal(unique(result$name), "TEST_SAMPLE")
})

test_that("read.ASAP.individual extracts numeric depth and breadth", {
  xml_path <- file.path(.fixture_dir, "sample_asap.xml")
  skip_if_not(file.exists(xml_path), "fixture XML not found")

  result <- suppressMessages(read.ASAP.individual(xml_path))

  expect_true(is.numeric(result$avg_depth))
  expect_true(is.numeric(result$breadth))
  expect_equal(result$avg_depth[1], 450)
  expect_equal(result$breadth[1],   0.98)
})

test_that("read.ASAP.snps.individual returns expected columns from fixture XML", {
  xml_path <- file.path(.fixture_dir, "sample_asap.xml")
  skip_if_not(file.exists(xml_path), "fixture XML not found")

  result <- suppressMessages(read.ASAP.snps.individual(xml_path))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)   # fixture has 2 SNP nodes

  expected_cols <- c("name", "assay_name", "snp_name", "snp_position",
                     "snp_reference", "snp_distribution", "snp_depth", "snp_proportion")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("read.ASAP.snps.individual parses snp_distribution string", {
  xml_path <- file.path(.fixture_dir, "sample_asap.xml")
  skip_if_not(file.exists(xml_path), "fixture XML not found")

  result <- suppressMessages(read.ASAP.snps.individual(xml_path))

  expect_true(all(grepl("=", result$snp_distribution)))
})
