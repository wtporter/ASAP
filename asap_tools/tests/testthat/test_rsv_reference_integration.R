# Integration tests against a REAL single-contig viral reference (RSVB) and,
# when available, real combined ASAP run data. Guards the current viral AA-calling
# path so the eukaryotic rework can't silently regress it.
#
# repo_root is derived from the functions dir (repo/asap_tools/asap_tools_functions).
.repo_root <- dirname(dirname(.functions_dir))
.rsvb_gb   <- file.path(.repo_root, "nextflow", "tests", "preparejson",
                        "RSVB_NC_001781.1.gb")
.rsv_rdata <- file.path(.repo_root, "nextflow", "tests", "test_output",
                        "RSV_ROI_Dev", "sample_reports", "rdata",
                        "RSV_ROI_Dev_ASAP_Data.Rdata")

test_that("extract_gene_table parses the real RSVB reference (11 CDS)", {
  skip_if_not(file.exists(.rsvb_gb), "RSVB reference not present")

  ref <- suppressWarnings(genbankr::readGenBank(.rsvb_gb))
  gt  <- extract_gene_table(ref)

  expect_equal(nrow(gt), 11L)
  expect_true(all(c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L")
                  %in% gt$gene))
  expect_false(any(is.na(gt$sequence)))
})

test_that("snps.to.amino calls a real RSVB in-gene SNP deterministically", {
  skip_if_not(file.exists(.rsvb_gb), "RSVB reference not present")

  # genome pos 781 sits in NS2 (CDS 626..1000); ref base T, mut A.
  # gene-relative (781-626)+1 = 156 -> synonymous at that codon.
  ref   <- suppressWarnings(genbankr::readGenBank(.rsvb_gb))
  token <- paste0(substr(as.character(ref@sequence), 781L, 781L), 781L, "A")
  snp_db <- data.frame(SNP = token, assay_name = "RSVB_NC_001781.1_NC_001781",
                       stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    snps.to.amino(snp_db, .rsvb_gb, cores = 1)))
  expect_equal(res$Gene, "NS2")
  expect_equal(res$SNP_Gene, "T156A")
  expect_equal(res$AA, "Synonymous")
})

test_that("genome.snp.to.gene.snp runs over real RSV SNP positions at scale", {
  skip_if_not(file.exists(.rsvb_gb),  "RSVB reference not present")
  skip_if_not(file.exists(.rsv_rdata), "RSV combined RData not present")

  load(.rsv_rdata)  # final_snps, final_array, final_asap
  ref <- suppressWarnings(genbankr::readGenBank(.rsvb_gb))

  rb <- final_snps[final_snps$assay_name == "RSVB_NC_001781.1_NC_001781", ]
  rb <- rb[!is.na(rb$snp_position) & !is.na(rb$snp_reference), ]
  # Build valid genome SNP tokens from real (reference, position) pairs.
  toks <- unique(paste0(rb$snp_reference, rb$snp_position, "A"))
  toks <- head(toks, 100L)  # cap for runtime
  snp_db <- data.frame(SNP = toks, assay_name = "RSVB_NC_001781.1_NC_001781",
                       stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    genome.snp.to.gene.snp(snp_db, .rsvb_gb, cores = 1)))

  # Every input SNP is accounted for, and RSV being gene-dense, most map to genes.
  expect_equal(nrow(res), length(toks))
  expect_true(all(c("SNP", "Gene", "SNP_Gene") %in% names(res)))
  expect_gt(sum(res$Gene != "Non-gene region"), 0L)
})
