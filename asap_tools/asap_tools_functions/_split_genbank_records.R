#' Split a (possibly multi-record) GenBank file into per-contig single-LOCUS files.
#'
#' genbankr::readGenBank() asserts exactly one LOCUS per file and errors on
#' multi-record input (bacterial/fungal assemblies, segmented viruses). This
#' splitter writes each record to its own single-LOCUS temp file so downstream
#' code can parse contigs one at a time.
#'
#' The LOCUS name captured here is the token ASAP uses to build assay_name
#' ("<filebase>_<LOCUS>", see asap/prepareJSONInput_nextflow.py), so it is the
#' key for matching SNPs to the right contig.
#'
#' @param gb_file  path to a GenBank file (1+ LOCUS records)
#' @param out_dir  directory for the per-record files (default: a fresh tempdir)
#' @return data.frame with columns: locus (chr), index (int, 1-based record
#'         order), path (chr, single-record file). One row per LOCUS record.
split_genbank_records <- function(gb_file, out_dir = tempfile("gbsplit_")) {
  lines <- readLines(gb_file, warn = FALSE)

  # Record boundaries: a record runs from a LOCUS line through its terminating
  # "//". Content outside those bounds (blank lines between records) is ignored.
  locus_idx <- grep("^LOCUS", lines)
  if (length(locus_idx) == 0L) {
    stop(sprintf("No LOCUS line found in GenBank file: %s", gb_file))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  file_base <- tools::file_path_sans_ext(basename(gb_file))
  end_idx   <- grep("^//", lines)

  recs <- vector("list", length(locus_idx))
  for (i in seq_along(locus_idx)) {
    start <- locus_idx[i]
    # this record ends at the first "//" at or after its LOCUS line
    stops <- end_idx[end_idx >= start]
    stop_line <- if (length(stops) > 0L) stops[1] else length(lines)

    locus_name <- strsplit(trimws(lines[start]), "\\s+")[[1]][2]
    safe_locus <- gsub("[^A-Za-z0-9._-]", "_", locus_name)
    out_path   <- file.path(out_dir, sprintf("%s__%s.gb", file_base, safe_locus))

    rec_lines <- lines[start:stop_line]

    # genbankr silently drops CDS features whose location has a partial marker
    # (`<`/`>`, e.g. `1..>2277` or `join(1..26,715..>979)`), yielding 0 CDS rows.
    # These markers appear only in feature location strings (always before a
    # digit), so removing them is safe and preserves reading frame (codon_start
    # is untouched).
    #
    # Before stripping, capture the genomic span (min..max of the location
    # integers) of each PARTIAL CDS. That span equals extract_gene_table()'s
    # start/end for the gene, so downstream code can mark exactly the truncated
    # gene(s) -- even on contigs that also carry a complete CDS (e.g. M2+M1 on
    # one influenza segment).
    partial_spans <- character(0)
    for (ci in grep("^     CDS ", rec_lines)) {          # CDS feature key line
      loc <- sub("^     CDS +", "", rec_lines[ci])
      j <- ci + 1L                                        # gather wrapped location lines
      while (j <= length(rec_lines) && grepl("^ {21}[^/]", rec_lines[j])) {
        loc <- paste0(loc, trimws(rec_lines[j])); j <- j + 1L
      }
      if (grepl("[<>]", loc)) {
        nums <- as.integer(regmatches(loc, gregexpr("[0-9]+", loc))[[1]])
        partial_spans <- c(partial_spans, paste(min(nums), max(nums)))
      }
    }
    rec_lines <- gsub("([<>])([0-9])", "\\2", rec_lines)

    writeLines(rec_lines, out_path)
    recs[[i]] <- data.frame(locus = locus_name, index = i, path = out_path,
                            stringsAsFactors = FALSE)
    recs[[i]]$partial_spans <- list(partial_spans)        # list-column, one row per record
  }

  do.call(rbind, recs)
}
