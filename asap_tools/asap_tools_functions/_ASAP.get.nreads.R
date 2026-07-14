library(dplyr)
library(foreach)
library(doParallel)

ASAP.get.nreads <- function(read.ASAP.df, num_cores = 1) {
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }

  n_reads_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run        <- read.ASAP.df$run[[i]]
    name       <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    n_reads    <- read.ASAP.df$n_reads[[i]]
    ref_pos    <- if ("ref_positions" %in% names(read.ASAP.df)) read.ASAP.df$ref_positions[[i]] else NA_character_

    n_reads  <- as.numeric(unlist(strsplit(n_reads, ",")))
    position <- as.numeric(unlist(strsplit(ref_pos, ",")))
    # Fall back to a contiguous 1-based index if ref_positions is absent/mismatched.
    if (length(position) != length(n_reads)) position <- seq_along(n_reads)

    data.frame(run, name, assay_name, position, n_reads)
  }

  n_reads_out$position <- as.numeric(as.character(n_reads_out$position))
  n_reads_out$n_reads  <- as.numeric(as.character(n_reads_out$n_reads))

  return(n_reads_out)
}
