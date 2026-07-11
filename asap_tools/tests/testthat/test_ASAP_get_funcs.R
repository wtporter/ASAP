# Tests for ASAP.get.depth, ASAP.get.proportions, ASAP.get.nreads,
# and ASAP.get.quality.discards.
#
# setup.R doesn't source these files, so we source them here.
skip_if_not_installed("foreach")
skip_if_not_installed("doParallel")

for (f in c("_ASAP.get.depth.R", "_ASAP.get.nreads.R",
            "_ASAP.get.proportions.R", "_ASAP.get.quality.discards.R")) {
  source(file.path(.functions_dir, f))
}

# Minimal read.ASAP.individual-shaped row for testing.
# ref_positions carries the genomic coordinate for each array entry (emitted 1:1 with the arrays,
# sparse under --prune-per-base). When it is length-matched to an array the extractor uses it for
# `position`; otherwise it falls back to a contiguous 1-based index.
make_asap_row <- function(run = "RUN1", name = "SAMPLE1", assay_name = "TB_rpoB",
                           depths = "100,200,300",
                           proportions = "0.95,0.98,0.97",
                           n_reads = "50,60,70",
                           quality_discards = "1,0,2",
                           ref_positions = "1,2,3") {
  data.frame(run, name, assay_name, depths, proportions,
             n_reads, quality_discards, ref_positions, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# ASAP.get.depth
# ---------------------------------------------------------------------------

test_that("ASAP.get.depth returns a data.frame", {
  result <- ASAP.get.depth(make_asap_row())
  expect_s3_class(result, "data.frame")
})

test_that("ASAP.get.depth has correct output columns", {
  result <- ASAP.get.depth(make_asap_row())
  expect_true(all(c("run", "name", "assay_name", "position", "depth") %in% names(result)))
})

test_that("ASAP.get.depth produces one row per comma-separated depth value", {
  result <- ASAP.get.depth(make_asap_row(depths = "100,200,300"))
  expect_equal(nrow(result), 3L)
})

test_that("ASAP.get.depth position is a 1-based integer sequence", {
  result <- ASAP.get.depth(make_asap_row(depths = "10,20,30"))
  expect_equal(result$position, c(1, 2, 3))
})

test_that("ASAP.get.depth depth values are numeric and match input", {
  result <- ASAP.get.depth(make_asap_row(depths = "10,20,30"))
  expect_true(is.numeric(result$depth))
  expect_equal(result$depth, c(10, 20, 30))
})

test_that("ASAP.get.depth handles a single-position amplicon", {
  result <- ASAP.get.depth(make_asap_row(depths = "450"))
  expect_equal(nrow(result), 1L)
  expect_equal(result$position, 1)
  expect_equal(result$depth, 450)
})

test_that("ASAP.get.depth stacks multiple amplicons without mixing positions", {
  df <- rbind(
    make_asap_row(assay_name = "amp1", depths = "10,20"),
    make_asap_row(assay_name = "amp2", depths = "30,40,50")
  )
  result <- ASAP.get.depth(df)
  expect_equal(nrow(result), 5L)
  expect_equal(sort(unique(result$assay_name)), c("amp1", "amp2"))
  # amp1 positions are 1-2, amp2 positions are 1-3
  expect_equal(result$position[result$assay_name == "amp1"], c(1, 2))
  expect_equal(result$position[result$assay_name == "amp2"], c(1, 2, 3))
})

test_that("ASAP.get.depth propagates run, name, assay_name to all rows", {
  result <- ASAP.get.depth(make_asap_row(run = "R1", name = "S1", assay_name = "GENE",
                                          depths = "5,10"))
  expect_true(all(result$run == "R1"))
  expect_true(all(result$name == "S1"))
  expect_true(all(result$assay_name == "GENE"))
})

test_that("ASAP.get.depth uses ref_positions as the genomic coordinate (sparse/--prune-per-base)", {
  # Pruned output: 3 covered positions scattered across the reference. position must be the
  # genomic coordinate from ref_positions, NOT a contiguous 1..3 index.
  result <- ASAP.get.depth(make_asap_row(depths = "100,150,200",
                                          ref_positions = "500,1200,4400"))
  expect_equal(result$position, c(500, 1200, 4400))
  expect_equal(result$depth, c(100, 150, 200))
})

test_that("ASAP.get.depth uses ref_positions with a non-1 genomic offset", {
  # Targeted assay whose reference starts at position 761101 (contiguous, but offset from 1).
  result <- ASAP.get.depth(make_asap_row(depths = "10,20,30",
                                          ref_positions = "761101,761102,761103"))
  expect_equal(result$position, c(761101, 761102, 761103))
})

test_that("ASAP.get.depth falls back to a 1-based index when ref_positions is absent", {
  row <- make_asap_row(depths = "10,20,30")
  row$ref_positions <- NULL   # e.g. Rdata predating --prune-per-base
  result <- ASAP.get.depth(row)
  expect_equal(result$position, c(1, 2, 3))
})

test_that("ASAP.get.depth falls back when ref_positions length mismatches the array", {
  result <- ASAP.get.depth(make_asap_row(depths = "10,20,30,40", ref_positions = "1,2,3"))
  expect_equal(result$position, c(1, 2, 3, 4))
})

test_that("ASAP.get.proportions uses ref_positions as the genomic coordinate", {
  result <- ASAP.get.proportions(make_asap_row(proportions = "0.9,0.8",
                                                ref_positions = "500,4400"))
  expect_equal(result$position, c(500, 4400))
})

# ---------------------------------------------------------------------------
# ASAP.get.proportions
# ---------------------------------------------------------------------------

test_that("ASAP.get.proportions returns a data.frame with correct columns", {
  result <- ASAP.get.proportions(make_asap_row())
  expect_s3_class(result, "data.frame")
  expect_true(all(c("run", "name", "assay_name", "position", "proportions") %in% names(result)))
})

test_that("ASAP.get.proportions produces one row per proportions value", {
  result <- ASAP.get.proportions(make_asap_row(proportions = "0.9,0.95,1.0"))
  expect_equal(nrow(result), 3L)
})

test_that("ASAP.get.proportions values are numeric and match input", {
  result <- ASAP.get.proportions(make_asap_row(proportions = "0.90,0.95,1.00"))
  expect_true(is.numeric(result$proportions))
  expect_equal(result$proportions, c(0.90, 0.95, 1.00))
})

test_that("ASAP.get.proportions position is 1-based", {
  result <- ASAP.get.proportions(make_asap_row(proportions = "0.9,0.8"))
  expect_equal(result$position, c(1, 2))
})

test_that("ASAP.get.proportions stacks multiple amplicons correctly", {
  df <- rbind(
    make_asap_row(assay_name = "a1", proportions = "0.9,0.8"),
    make_asap_row(assay_name = "a2", proportions = "1.0")
  )
  result <- ASAP.get.proportions(df)
  expect_equal(nrow(result), 3L)
})

# ---------------------------------------------------------------------------
# ASAP.get.nreads
# ---------------------------------------------------------------------------

test_that("ASAP.get.nreads returns a data.frame with correct columns", {
  result <- ASAP.get.nreads(make_asap_row())
  expect_s3_class(result, "data.frame")
  expect_true(all(c("run", "name", "assay_name", "position", "n_reads") %in% names(result)))
})

test_that("ASAP.get.nreads produces one row per n_reads value", {
  result <- ASAP.get.nreads(make_asap_row(n_reads = "50,60,70"))
  expect_equal(nrow(result), 3L)
})

test_that("ASAP.get.nreads n_reads values are numeric and match input", {
  result <- ASAP.get.nreads(make_asap_row(n_reads = "10,20,30"))
  expect_true(is.numeric(result$n_reads))
  expect_equal(result$n_reads, c(10, 20, 30))
})

test_that("ASAP.get.nreads position is 1-based", {
  result <- ASAP.get.nreads(make_asap_row(n_reads = "5,10"))
  expect_equal(result$position, c(1, 2))
})

test_that("ASAP.get.nreads stacks multiple amplicons correctly", {
  df <- rbind(
    make_asap_row(assay_name = "amp1", n_reads = "10,20"),
    make_asap_row(assay_name = "amp2", n_reads = "30,40,50")
  )
  result <- ASAP.get.nreads(df)
  expect_equal(nrow(result), 5L)
  expect_equal(sort(unique(result$assay_name)), c("amp1", "amp2"))
})

# ---------------------------------------------------------------------------
# ASAP.get.quality.discards
# ---------------------------------------------------------------------------

test_that("ASAP.get.quality.discards returns a data.frame with correct columns", {
  result <- ASAP.get.quality.discards(make_asap_row())
  expect_s3_class(result, "data.frame")
  expect_true(all(c("run", "name", "assay_name", "position", "quality_discards") %in% names(result)))
})

test_that("ASAP.get.quality.discards produces one row per quality_discards value", {
  result <- ASAP.get.quality.discards(make_asap_row(quality_discards = "1,0,2"))
  expect_equal(nrow(result), 3L)
})

test_that("ASAP.get.quality.discards values are numeric and match input", {
  result <- ASAP.get.quality.discards(make_asap_row(quality_discards = "1,0,2"))
  expect_true(is.numeric(result$quality_discards))
  expect_equal(result$quality_discards, c(1, 0, 2))
})

test_that("ASAP.get.quality.discards position is 1-based", {
  result <- ASAP.get.quality.discards(make_asap_row(quality_discards = "0,1,2,3"))
  expect_equal(result$position, c(1, 2, 3, 4))
})

test_that("ASAP.get.quality.discards handles all-zero discards (clean amplicon)", {
  result <- ASAP.get.quality.discards(make_asap_row(quality_discards = "0,0,0"))
  expect_equal(nrow(result), 3L)
  expect_equal(result$quality_discards, c(0, 0, 0))
})

test_that("ASAP.get.quality.discards stacks multiple amplicons correctly", {
  df <- rbind(
    make_asap_row(assay_name = "amp1", quality_discards = "0,1"),
    make_asap_row(assay_name = "amp2", quality_discards = "2,3,4")
  )
  result <- ASAP.get.quality.discards(df)
  expect_equal(nrow(result), 5L)
})
