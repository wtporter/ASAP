library(dplyr)
library(foreach)
library(doParallel)

ASAP.get.depth <- function(read.ASAP.df, num_cores = 1) {
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }

  depth_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run        <- read.ASAP.df$run[[i]]
    name       <- read.ASAP.df$name[[i]]
    assay_name <- read.ASAP.df$assay_name[[i]]
    depth      <- read.ASAP.df$depths[[i]]

    depth    <- as.numeric(unlist(strsplit(depth, ",")))
    position <- seq_along(depth)

    data.frame(run, name, assay_name, position, depth)
  }

  depth_out$position <- as.numeric(as.character(depth_out$position))
  depth_out$depth    <- as.numeric(as.character(depth_out$depth))

  return(depth_out)
}
