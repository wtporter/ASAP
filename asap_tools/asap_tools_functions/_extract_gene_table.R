library(genbankr)

#' Build a per-CDS gene table from a genbankr reference, defensively.
#'
#' Always derived from `reference@cds` (present for every GenBank file) rather
#' than `reference@genes` (built from `gene`-type features, which some viral
#' RefSeq records omit entirely -- e.g. CDS-only annotations with no `/gene`
#' or `/locus_tag` qualifier, only `/product`). Guarantees `gene`, `product`,
#' `translation`, and `sequence` columns exist regardless of what the source
#' file actually annotated.
extract_gene_table <- function(reference) {
  df <- data.frame(reference@cds)

  # Flatten any CharacterList columns (e.g. translation) to plain strings.
  df <- as.data.frame(lapply(df, function(x) if (is.list(x)) sapply(x, paste, collapse = ";") else x))

  # Label fallback chain: gene -> locus_tag -> product -> generated placeholder.
  # gene is the short, conventional symbol (e.g. "NS1") when annotated; CDS-only
  # files have neither a `gene` nor a `locus_tag` column at all -- only
  # `product`, which is often already a real gene name (e.g. "E1A").
  for (col in c("gene", "locus_tag", "product")) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  df$gene <- with(df, ifelse(!is.na(gene), gene,
                       ifelse(!is.na(locus_tag), locus_tag,
                       ifelse(!is.na(product), product,
                              paste0("CDS_", seq_len(nrow(df)))))))

  if (!"translation" %in% names(df)) df$translation <- NA_character_

  df$sequence <- "No Seq"
  for (i in seq_len(nrow(df))) {
    df$sequence[i] <- substr(as.character(reference@sequence), df$start[i], df$end[i])
  }

  df
}
