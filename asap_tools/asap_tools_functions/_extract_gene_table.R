library(genbankr)

#' Build a per-transcript CDS table from a genbankr reference, defensively.
#'
#' Derived from `reference@cds`. genbankr emits ONE ROW PER EXON for a spliced
#' (join/order) CDS, all sharing a `transcript_id`; this collapses each
#' transcript to a single record whose `sequence` is the SPLICED coding region
#' (introns excised, exons concatenated in ascending genomic / + strand order)
#' and whose `exons` list-column carries the exon coordinates needed to map a
#' genomic SNP position onto the spliced CDS. Single-exon (viral) CDS collapse
#' to one row identical to the pre-splicing behavior.
#'
#' Label fallback chain: gene -> locus_tag -> product -> generated placeholder,
#' so `gene`, `product`, `sequence`, `strand`, `start`, `end`, `codon_start`,
#' `translation`, and `exons` always exist regardless of what the source
#' annotated (CDS-only viral records omit `gene`/`locus_tag`).
#'
#' @return data.frame, one row per transcript, with columns gene, product,
#'   strand ("+"/"-"), start (min exon), end (max exon), codon_start,
#'   translation, sequence (spliced), and exons (list-column of data.frame
#'   with start/end, ascending).
extract_gene_table <- function(reference) {
  df <- data.frame(reference@cds)

  # Flatten any CharacterList columns (e.g. translation) to plain strings.
  df <- as.data.frame(lapply(df, function(x) if (is.list(x)) sapply(x, paste, collapse = ";") else x),
                      stringsAsFactors = FALSE)

  for (col in c("gene", "locus_tag", "product")) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  if (!"translation" %in% names(df)) df$translation <- NA_character_
  if (!"codon_start" %in% names(df)) df$codon_start <- 1L

  # Transcript grouping key: transcript_id (present even for CDS-only records),
  # then gene_id, then a unique per-row id so distinct un-keyed CDS never merge.
  key <- rep(NA_character_, nrow(df))
  if ("transcript_id" %in% names(df)) key <- as.character(df$transcript_id)
  if ("gene_id" %in% names(df)) { na <- is.na(key); key[na] <- as.character(df$gene_id)[na] }
  na <- is.na(key); key[na] <- paste0(".ungrouped_", seq_len(nrow(df)))[na]

  full_seq <- as.character(reference@sequence)
  groups   <- unique(key)  # first-appearance order

  # first non-NA / non-empty value of a column within a transcript's exon rows
  pick <- function(sub, col) {
    if (!col %in% names(sub)) return(NA_character_)
    v <- sub[[col]][!is.na(sub[[col]]) & nzchar(as.character(sub[[col]]))]
    if (length(v)) as.character(v[1]) else NA_character_
  }

  rows <- lapply(groups, function(g) {
    sub   <- df[key == g, , drop = FALSE]
    sub   <- sub[order(sub$start), , drop = FALSE]
    exons <- data.frame(start = as.integer(sub$start),
                        end   = as.integer(sub$end),
                        stringsAsFactors = FALSE)

    spliced <- paste0(mapply(function(s, e) substr(full_seq, s, e),
                             exons$start, exons$end), collapse = "")

    gene <- pick(sub, "gene")
    if (is.na(gene)) gene <- pick(sub, "locus_tag")
    if (is.na(gene)) gene <- pick(sub, "product")
    if (is.na(gene)) gene <- paste0("CDS_", g)

    out <- data.frame(
      gene        = gene,
      product     = pick(sub, "product"),
      strand      = as.character(sub$strand[1]),
      start       = min(exons$start),
      end         = max(exons$end),
      codon_start = suppressWarnings(as.integer(sub$codon_start[1])),
      translation = pick(sub, "translation"),
      sequence    = spliced,
      stringsAsFactors = FALSE
    )
    out$exons <- list(exons)   # per-transcript exon coordinates
    out
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Map a genomic (1-based, + strand) position onto a spliced CDS coordinate.
#'
#' Excises introns: returns the 1-based offset of `pos` within the concatenated
#' exon (+ strand) sequence, or NA if `pos` falls in an intron or outside the
#' transcript. `exons` is a data.frame(start, end) in ascending genomic order
#' (as produced by extract_gene_table). For a single-exon gene this reduces to
#' (pos - start) + 1, matching the pre-splicing behavior.
genomic_to_spliced_pos <- function(pos, exons) {
  w <- which(exons$start <= pos & exons$end >= pos)
  if (length(w) == 0L) return(NA_integer_)
  w <- w[1]
  prior <- if (w > 1L)
    sum(exons$end[seq_len(w - 1L)] - exons$start[seq_len(w - 1L)] + 1L) else 0L
  as.integer(prior + (pos - exons$start[w] + 1L))
}
