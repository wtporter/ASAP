library(dplyr)
library(foreach)
library(doParallel)

ASAP.get.proportions <- function(read.ASAP.df, num_cores = 1) {
  if (num_cores > 1) {
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }

  proportions_out <- foreach(i = 1:nrow(read.ASAP.df), .combine = bind_rows) %dopar% {
    run         <- read.ASAP.df$run[[i]]
    name        <- read.ASAP.df$name[[i]]
    assay_name  <- read.ASAP.df$assay_name[[i]]
    proportions <- read.ASAP.df$proportions[[i]]

    proportions <- as.numeric(unlist(strsplit(proportions, ",")))
    position    <- seq_along(proportions)

    data.frame(run, name, assay_name, position, proportions)
  }

  proportions_out$position    <- as.numeric(as.character(proportions_out$position))
  proportions_out$proportions <- as.numeric(as.character(proportions_out$proportions))

  return(proportions_out)
}
