# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)
library(ggplot2)

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
  if (solver1$name == "rk45") {
    rel_solvers <- rk45_list(tols = confs)
  } else if (solver1$name == "rk4") {
    rel_solvers <- rk4_list(num_steps = confs)
  } else if (solver1$name == "midpoint") {
    rel_solvers <- midpoint_list(num_steps = confs)
  } else {
    stop("unknown solver")
  }

  # Run reliability check
  reliab <- fit$reliability(
    solvers = rel_solvers, force = TRUE, savedir = dir
  )

  # Return list
  list(fits = fits, reliab = reliab, confs = confs, res = res)
}

# Run
r1 <- reliability("results_rk45", 3)
r2 <- reliability("results_rk4", 3)
r3 <- reliability("results_midpoint", 3)

# Plots
#plt1 <- plot_metrics(reliab, num_steps = ns_rel)
#plt2 <- plot_time_comparison_ns(fits, reliab, idx)
#plt3 <- plot_time_comparison_ns(fits, reliab, idx, TRUE)
#diags <- get_diags_df(fits) # rhat and reff
