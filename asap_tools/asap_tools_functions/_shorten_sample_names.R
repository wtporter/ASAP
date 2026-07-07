#' Shorten sample names for plot labels by stripping boilerplate.
#'
#' Two steps, applied across the full set of unique names so every sample in
#' a cohort gets a consistently-derived label. Callers can pass a raw column
#' with repeated names (e.g. one row per sample per amplicon) -- duplicates
#' are collapsed internally before computing the shortening, and every
#' occurrence of a given name maps back to the same shortened label.
#'  1. Strip the Illumina sample-index/lane suffix (`_S<digits>_L...`). The
#'     match requires `_L` right after the digits so an incidental "S<digits>"
#'     elsewhere in a sample name isn't mistaken for it.
#'  2. Strip the longest common prefix left across all names, snapped back to
#'     the nearest `_` boundary so a shared token isn't cut mid-word (e.g.
#'     "ADV_HAdV41_CoF_Rio..." / "ADV_HAdV41_CoT_Tempe..." keeps "CoF"/"CoT"
#'     intact rather than trimming to "F_Rio..." / "T_Tempe...").
#'
#' If either step would leave an empty or non-unique label for a given
#' sample, that sample's original (untrimmed) name is kept instead. As a
#' final guarantee, the returned labels are always unique and non-blank --
#' if two *distinct* input names still collide even after falling back to
#' their originals (i.e. the originals themselves were identical), a numeric
#' suffix is appended so labels never silently duplicate.
shorten_sample_names <- function(names) {
  unique_names <- unique(names)

  stripped <- sub("_S[0-9]+_L.*$", "", unique_names)

  common_prefix_len <- 0L
  if (length(stripped) > 1) {
    min_len <- min(nchar(stripped))
    if (min_len > 0) {
      for (i in seq_len(min_len)) {
        chars_at_i <- substr(stripped, i, i)
        if (length(unique(chars_at_i)) == 1) {
          common_prefix_len <- i
        } else {
          break
        }
      }
    }
  }
  # Snap the cut point back to the last "_" at or before common_prefix_len,
  # so we never split a shared token in half.
  if (common_prefix_len > 0) {
    prefix <- substr(stripped[1], 1, common_prefix_len)
    last_underscore <- max(gregexpr("_", prefix, fixed = TRUE)[[1]])
    if (last_underscore > 0) common_prefix_len <- last_underscore
  }
  shortened <- substr(stripped, common_prefix_len + 1, nchar(stripped))

  bad <- shortened == "" | duplicated(shortened) | duplicated(shortened, fromLast = TRUE)
  shortened[bad] <- unique_names[bad]

  # Final guarantee: unique, non-blank labels no matter what.
  still_bad <- shortened == "" | duplicated(shortened) | duplicated(shortened, fromLast = TRUE)
  if (any(still_bad)) {
    shortened[still_bad] <- make.unique(shortened[still_bad], sep = "_")
  }

  setNames(shortened, unique_names)[names]
}
