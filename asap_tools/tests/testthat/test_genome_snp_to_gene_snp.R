# Tests for genome.snp.to.gene.snp() — maps a genome-coordinate SNP token
# ("<ref><pos><mut>") onto the gene it falls in and returns the gene-relative
# SNP string. Previously had no direct test coverage.
#
# Uses the single-exon NS1 fixture (tiny_with_locus_tag.gb): CDS 99..518, +strand.

snp_db_in_gene <- function(pos = 150L, mut = "A") {
  ref <- suppressWarnings(genbankr::readGenBank(
    file.path(.fixture_dir, "tiny_with_locus_tag.gb")))
  refbase <- substr(as.character(ref@sequence), pos, pos)
  data.frame(SNP        = paste0(refbase, pos, mut),
             assay_name = "TESTSEQ1_NS1",
             stringsAsFactors = FALSE)
}

test_that("genome.snp.to.gene.snp maps an in-gene SNP to gene coordinates", {
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  res <- suppressWarnings(suppressMessages(
    genome.snp.to.gene.snp(snp_db_in_gene(150L, "A"), gb, cores = 1)))

  expect_equal(nrow(res), 1L)
  expect_equal(res$Gene, "NS1")
  # genome pos 150 in CDS starting at 99 -> gene-relative (150-99)+1 = 52
  expect_equal(res$SNP_Gene, "G52A")
})

test_that("genome.snp.to.gene.snp returns the documented output columns", {
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  res <- suppressWarnings(suppressMessages(
    genome.snp.to.gene.snp(snp_db_in_gene(150L, "A"), gb, cores = 1)))
  expect_true(all(c("SNP", "Gene", "SNP_Gene") %in% names(res)))
})

test_that("genome.snp.to.gene.snp labels out-of-gene SNPs as Non-gene region (mixed batch)", {
  # Batch with one in-gene (150) and one intergenic (20, before CDS start 99).
  # The intergenic SNP must come back as 'Non-gene region', not error.
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  ref <- suppressWarnings(genbankr::readGenBank(gb))
  mk  <- function(p) paste0(substr(as.character(ref@sequence), p, p), p, "A")
  snp_db <- data.frame(SNP = c(mk(150L), mk(20L)),
                       assay_name = "TESTSEQ1_NS1", stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    genome.snp.to.gene.snp(snp_db, gb, cores = 1)))
  expect_equal(nrow(res), 2L)
  expect_equal(res$Gene[res$SNP == mk(150L)], "NS1")
  expect_equal(res$Gene[res$SNP == mk(20L)],  "Non-gene region")
})

test_that("genome.snp.to.gene.snp handles an all-intergenic batch without erroring", {
  # Regression: when NO SNP in the batch maps to any gene, the internal foreach
  # returns an empty frame with no columns; the final full_join() used to abort
  # with "`by` must be supplied ...". A batch of purely intronic SNPs against a
  # spliced reference hits exactly this. Desired: 'Non-gene region' rows.
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  ref <- suppressWarnings(genbankr::readGenBank(gb))
  snp <- paste0(substr(as.character(ref@sequence), 20L, 20L), 20L, "A")
  snp_db <- data.frame(SNP = snp, assay_name = "TESTSEQ1_NS1",
                       stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    genome.snp.to.gene.snp(snp_db, gb, cores = 1)))
  expect_equal(nrow(res), 1L)
  expect_equal(res$Gene, "Non-gene region")
})
