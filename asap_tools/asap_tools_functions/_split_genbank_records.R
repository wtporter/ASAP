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

    writeLines(lines[start:stop_line], out_path)
    recs[[i]] <- data.frame(locus = locus_name, index = i, path = out_path,
                            stringsAsFactors = FALSE)
  }

  do.call(rbind, recs)
}
