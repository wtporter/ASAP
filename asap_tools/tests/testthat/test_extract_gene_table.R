test_that("extract_gene_table uses gene/locus_tag when annotated", {
  ref <- suppressWarnings(genbankr::readGenBank(
    file.path(.fixture_dir, "tiny_with_locus_tag.gb")))
  result <- extract_gene_table(ref)

  expect_equal(nrow(result), 1L)
  expect_equal(result$gene[1], "NS1")
  expect_equal(result$start[1], 99L)
  expect_equal(result$end[1], 518L)
  expect_false(is.na(result$sequence[1]))
  expect_false(any(is.na(result$gene)))
})

test_that("extract_gene_table falls back to product when gene/locus_tag are absent", {
  # Mirrors the real bug: some GenBank records (e.g. CDS-only viral annotations)
  # have no `gene`-type features at all and no /locus_tag qualifier on CDS,
  # which used to make `select(locus_tag, ...)` hard-error.
  ref <- suppressWarnings(genbankr::readGenBank(
    file.path(.fixture_dir, "tiny_cds_only.gb")))
  result <- extract_gene_table(ref)

  expect_equal(nrow(result), 1L)
  expect_equal(result$gene[1], "nonstructural protein 1")
  expect_equal(result$product[1], "nonstructural protein 1")
  expect_false(is.na(result$sequence[1]))
  expect_false(any(is.na(result$gene)))
})
