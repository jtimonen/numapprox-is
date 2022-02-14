# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Define function
reliability <- function(dir, idx) {
  res <- readRDS(file.path(dir, "sampling.rds"))
  fits <- res$fits
  solver1 <- fits$solvers[[1]]
  if (is(solver1, "AdaptiveOdeSolver")) {
    confs <- get_tol_vec(fits$solvers)
  } else {
    confs <- get_num_steps_vec(fits$solvers)
  }

  # Load fit
  inds_rel <- (1 + idx):length(confs)
  fit <- load_fit(file = fits$files[idx])
  confs_rel <- confs[inds_rel]
  if (solver1$name == "bdf") {
    rel_solvers <- bdf_list(tols = confs_rel)
  } else {
    stop("unknown solver")
  }

  # Run reliability check
  bn <- paste0("odegq_idx", idx)
  reliab <- fit$reliability(
    solvers = rel_solvers, force = TRUE, savedir = dir, basename = bn
  )

  # Return list
  list(
    reliab = reliab,
    confs = confs,
    confs_rel = confs_rel,
    res = res,
    idx = idx
  )
}

# Run
res_dir <- c("results_bdf")
start_inds <- c(1, 2, 3, 4)
outputs <- list()
j <- 0
for (index in start_inds) {
  j <- j + 1
  outputs[[j]] <- reliability(res_dir, index)
}
fp <- file.path(res_dir, "reliability.rds")
out <- list(outputs = outputs, start_inds = start_inds)
saveRDS(out, file = fp)
