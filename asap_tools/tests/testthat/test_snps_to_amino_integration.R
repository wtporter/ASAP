# End-to-end tests for the snps.to.amino() ORCHESTRATOR (not just the pure
# .translate_* helpers covered in test_snps_to_amino.R). Exercises the SNP-token
# tokenizer, the SNP/insertion/deletion classification, the foreach dispatch,
# and the final join/labelling — against the single-exon NS1 fixture.

amino_on_ns1 <- function(snp_token) {
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  snp_db <- data.frame(SNP = snp_token, assay_name = "TESTSEQ1_NS1",
                       stringsAsFactors = FALSE)
  suppressWarnings(suppressMessages(snps.to.amino(snp_db, gb, cores = 1)))
}

test_that("snps.to.amino calls a missense substitution end-to-end", {
  # genome pos 150 (ref G) -> gene-relative 52 -> codon 18; G->A gives D18N
  res <- amino_on_ns1("G150A")
  expect_equal(nrow(res), 1L)
  expect_equal(res$Gene, "NS1")
  expect_equal(res$SNP_Gene, "G52A")
  expect_equal(res$AA, "NS1:D18N")
  expect_equal(res$Product, "nonstructural protein 1")
})

test_that("snps.to.amino returns the documented output columns", {
  res <- amino_on_ns1("G150A")
  expect_true(all(c("SNP", "SNP_Gene", "AA", "Gene", "Product",
                    "Theoretical_Reference") %in% names(res)))
})

test_that("snps.to.amino labels an intergenic substitution Non-coding SNP (mixed batch)", {
  # Batch of one in-gene (150) + one intergenic (20, upstream of CDS start 99).
  gb  <- file.path(.fixture_dir, "tiny_with_locus_tag.gb")
  ref <- suppressWarnings(genbankr::readGenBank(gb))
  mk  <- function(p) paste0(substr(as.character(ref@sequence), p, p), p, "A")
  snp_db <- data.frame(SNP = c(mk(150L), mk(20L)),
                       assay_name = "TESTSEQ1_NS1", stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(snps.to.amino(snp_db, gb, cores = 1)))
  expect_equal(res$AA[res$SNP == mk(20L)], "Non-coding SNP")
  expect_equal(res$Gene[res$SNP == mk(150L)], "NS1")
})

test_that("snps.to.amino handles an all-intergenic batch without erroring", {
  # Regression: identical to the genome.snp.to.gene.snp crash — when NO SNP maps
  # to a gene, the empty foreach result made the final full_join() abort with
  # "`by` must be supplied ...". Desired: a single 'Non-coding SNP' row.
  # All-intronic batches against spliced refs hit this.
  ref <- suppressWarnings(genbankr::readGenBank(
    file.path(.fixture_dir, "tiny_with_locus_tag.gb")))
  token <- paste0(substr(as.character(ref@sequence), 20L, 20L), 20L, "A")
  res <- amino_on_ns1(token)
  expect_equal(nrow(res), 1L)
  expect_equal(res$AA, "Non-coding SNP")
})

test_that("snps.to.amino classifies an in-frame deletion", {
  # 3-bp deletion inside the CDS: positions 150-152 deleted (mut '_')
  # Deletion tokens are '<ref><pos>_' for each deleted base, pipe-joined.
  ref <- suppressWarnings(genbankr::readGenBank(
    file.path(.fixture_dir, "tiny_with_locus_tag.gb")))
  bases <- vapply(150:152, function(p)
    substr(as.character(ref@sequence), p, p), character(1))
  token <- paste(paste0(bases, 150:152, "_"), collapse = "|")
  res <- amino_on_ns1(token)
  expect_equal(nrow(res), 1L)
  expect_equal(res$Gene, "NS1")
  expect_false(is.na(res$AA))
  # in-frame (3 bp) deletion should not be flagged as a frameshift
  expect_false(res$AA %in% c("Deletion Not In-frame", "Non-coding SNP"))
})
