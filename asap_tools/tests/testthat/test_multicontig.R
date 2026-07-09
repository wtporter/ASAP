# Tests for multi-contig (multi-LOCUS) GenBank handling.
#
# Fixture tiny_multicontig.gb concatenates two single-LOCUS records
# (TESTSEQ1 with the NS1 CDS, TESTSEQ2). It reflects the real eukaryotic case:
# ASAP names assays "<filebase>_<LOCUS>" per contig.

MULTICONTIG_GB <- function() file.path(.fixture_dir, "tiny_multicontig.gb")

test_that("genbankr cannot read a multi-record file (motivates record splitting)", {
  # Documents the hard constraint driving the planned splitter: readGenBank
  # asserts a single LOCUS per file and errors on multi-record input. This is
  # why per-contig temp files must be produced before parsing.
  expect_error(
    suppressWarnings(genbankr::readGenBank(MULTICONTIG_GB())),
    regexp = "LOCUS"
  )
})

test_that("split_genbank_records yields one single-LOCUS record per contig", {
  # split a multi-record GenBank into per-contig entries carrying the LOCUS name
  # (the token ASAP uses to build assay_name) and a parseable single-record path.
  recs <- split_genbank_records(MULTICONTIG_GB())

  expect_equal(nrow(recs), 2L)
  expect_setequal(recs$locus, c("TESTSEQ1", "TESTSEQ2"))
  # each split record must itself be readable by genbankr (no multi-LOCUS error)
  for (p in recs$path) {
    expect_error(suppressWarnings(genbankr::readGenBank(p)), NA)
  }
})

test_that("split_genbank_records passes a single-record file through unchanged", {
  recs <- split_genbank_records(
    file.path(.fixture_dir, "tiny_with_locus_tag.gb"))
  expect_equal(nrow(recs), 1L)
  expect_equal(recs$locus, "TESTSEQ1")
})

test_that("per-contig SNP scoping does not map across contigs", {
  # The pipeline splits the multi-contig file and translates each contig's SNPs
  # against ONLY that contig's record. TESTSEQ1 carries NS1 (CDS 99..518);
  # TESTSEQ2 carries its own CDS at the same coordinates. A SNP scoped to
  # TESTSEQ2 must map to TESTSEQ2's gene, never NS1.
  recs <- split_genbank_records(MULTICONTIG_GB())
  p2   <- recs$path[recs$locus == "TESTSEQ2"]
  ref2 <- suppressWarnings(genbankr::readGenBank(p2))
  token <- paste0(substr(as.character(ref2@sequence), 150L, 150L), 150L, "A")
  snp_db <- data.frame(SNP = token, assay_name = "TESTSEQ2_TESTSEQ2",
                       stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    snps.to.amino(snp_db, p2, cores = 1)))
  expect_false("NS1" %in% res$Gene)                 # TESTSEQ1's gene never surfaces
  expect_true(any(grepl("nonstructural", res$Gene))) # maps to TESTSEQ2's own CDS
})
