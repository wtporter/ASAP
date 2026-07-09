# TARGET-SPEC tests for spliced (multi-exon / intron-containing) CDS handling.
#
# These encode the CORRECT behavior that eukaryotic AA calling requires but that
# the current single-interval code does not yet implement. They are skipped so
# the suite stays green; un-skip them as the exon-aware rework lands. See
# project_eukaryotic_aa_calling in the session memory for the full analysis.
#
# Fixture tiny_spliced_cds.gb: CDS join(1..8,21..36) on +strand, one intron at
# genomic 9..20. Spliced CDS = ATGGCAAA + ACCCGGGTTTTGGTAA = ATGGCAAAACCCGGGTTTTGGTAA
# which translates to MAKPGFW (verified with Biostrings::translate).

SPLICED_GB <- function() file.path(.fixture_dir, "tiny_spliced_cds.gb")

test_that("extract_gene_table splices a multi-exon CDS into one gene", {
  ref <- suppressWarnings(genbankr::readGenBank(SPLICED_GB()))
  gt  <- extract_gene_table(ref)

  # One spliced gene, not two exon fragments
  expect_equal(nrow(gt), 1L)
  expect_equal(gt$gene[1], "SPLG")
  # Intron (genomic 9..20) excised: spliced CDS is exon1 + exon2, 24 bp
  expect_equal(gt$sequence[1], "ATGGCAAAACCCGGGTTTTGGTAA")
})

test_that("snps.to.amino intron-corrects an exonic SNP position", {
  # Genomic 25 (ref G) sits in exon 2. Excising the intron, it is spliced base 13
  # = codon 5 position 1; G->T turns codon 5 (GGG=Gly) into TGG=Trp -> G5W.
  # The current per-exon math mis-numbers this because it ignores the intron.
  snp_db <- data.frame(SNP = "G25T", assay_name = "SPLICED1_SPLG",
                       stringsAsFactors = FALSE)
  res <- suppressWarnings(suppressMessages(
    snps.to.amino(snp_db, SPLICED_GB(), cores = 1)))

  expect_equal(res$Gene, "SPLG")
  expect_equal(res$AA,   "SPLG:G5W")
})

test_that("snps.to.amino labels an intronic SNP as non-coding", {
  # Genomic 15 falls inside the intron (9..20): not part of the spliced CDS.
  snp_db <- data.frame(SNP = "A15G", assay_name = "SPLICED1_SPLG",
                       stringsAsFactors = FALSE)
  res <- suppressWarnings(suppressMessages(
    snps.to.amino(snp_db, SPLICED_GB(), cores = 1)))

  expect_match(res$AA, "[Nn]on-coding")
})
