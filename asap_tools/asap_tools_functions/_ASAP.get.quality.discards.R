library(dplyr)
library(foreach)
library(doParallel)

ASAP.get.quality.discards <- function(read.ASAP.df, num_cores = 1) {
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }

  quality_discards_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run              <- read.ASAP.df$run[[i]]
    name             <- read.ASAP.df$name[[i]]
    assay_name       <- read.ASAP.df$assay_name[[i]]
    quality_discards <- read.ASAP.df$quality_discards[[i]]
    ref_pos          <- if ("ref_positions" %in% names(read.ASAP.df)) read.ASAP.df$ref_positions[[i]] else NA_character_

    quality_discards <- as.numeric(unlist(strsplit(quality_discards, ",")))
    position         <- as.numeric(unlist(strsplit(ref_pos, ",")))
    # Fall back to a contiguous 1-based index if ref_positions is absent/mismatched.
    if (length(position) != length(quality_discards)) position <- seq_along(quality_discards)

    data.frame(run, name, assay_name, position, quality_discards)
  }

  quality_discards_out$position         <- as.numeric(as.character(quality_discards_out$position))
  quality_discards_out$quality_discards <- as.numeric(as.character(quality_discards_out$quality_discards))

  return(quality_discards_out)
}
