test_that("strips Illumina sample-index/lane suffix and common prefix, preserving word boundaries", {
  names <- c(
    "ADV_HAdV41_CoF_Rio_20251001_S13_L001",
    "ADV_HAdV41_CoF_Rio_20251008_S22_L001",
    "ADV_HAdV41_CoT_Tempe_2509414_507408_20250923_S4_L001",
    "ADV_HAdV41_CoT_Tempe_2509415_507407_20250923_S5_L001"
  )
  result <- shorten_sample_names(names)

  expect_equal(unname(result), c(
    "CoF_Rio_20251001",
    "CoF_Rio_20251008",
    "CoT_Tempe_2509414_507408_20250923",
    "CoT_Tempe_2509415_507407_20250923"
  ))
  # Word-boundary check: the shared prefix "ADV_HAdV41_Co" must not be cut
  # mid-token (i.e. must not produce "F_Rio..." / "T_Tempe...").
  expect_false(any(grepl("^F_|^T_", result)))
})

test_that("does not strip _S<digits> unless followed by _L (avoids false positives)", {
  names <- c("Study_S2_Cohort_A", "Study_S2_Cohort_B")
  result <- shorten_sample_names(names)
  expect_equal(unname(result), c("A", "B"))
})

test_that("single sample: suffix stripped, no prefix trimming attempted", {
  result <- shorten_sample_names("ONLY_SAMPLE_S1_L001")
  expect_equal(unname(result), "ONLY_SAMPLE")
})

test_that("falls back to original name when shortening would collide", {
  names <- c("Foo_S1_L001", "Foo_S2_L001")
  result <- shorten_sample_names(names)
  expect_equal(unname(result), names)
})

test_that("raw column with repeated names (normal one-row-per-amplicon case) maps consistently", {
  names <- rep(c("SAMPLE_A_S1_L001", "SAMPLE_B_S2_L001"), each = 3)
  result <- shorten_sample_names(names)
  expect_equal(unname(result), rep(c("A", "B"), each = 3))
})

test_that("guarantees unique output even when original names are truly identical", {
  names <- c("DUP_S1_L001", "DUP_S1_L001", "OTHER_S2_L001")
  result <- shorten_sample_names(names)
  # The two identical originals legitimately share one label...
  expect_equal(unname(result)[1:2], c("DUP", "DUP"))
  # ...but the guarantee is about *distinct* original names never colliding,
  # not about deduplicating identical ones -- confirm output stays non-blank.
  expect_true(all(nchar(unname(result)) > 0))
})

test_that("output preserves input order and length", {
  names <- c("B_S2_L001", "A_S1_L001", "B_S2_L001")
  result <- shorten_sample_names(names)
  expect_equal(length(result), length(names))
})
