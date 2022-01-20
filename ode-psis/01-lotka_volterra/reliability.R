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
    rel_solvers <- rk45_list(tols = confs_rel)
  } else if (solver1$name == "rk4") {
    rel_solvers <- rk4_list(num_steps = confs_rel)
  } else if (solver1$name == "midpoint") {
    rel_solvers <- midpoint_list(num_steps = confs_rel)
  } else {
    stop("unknown solver")
  }

  # Run reliability check
  reliab <- fit$reliability(
    solvers = rel_solvers, force = TRUE, savedir = dir
  )

  # Return list
  list(
    fits = fits,
    reliab = reliab,
    confs = confs,
    confs_rel = confs_rel,
    res = res,
    idx = idx
  )
}

# Run
idx_1 <- 3 # 1 and 2 fail
idx_2 <- 1
idx_3 <- 1
r1 <- reliability("results_rk45", idx_1)
r2 <- reliability("results_rk4", idx_2)
r3 <- reliability("results_midpoint", idx_3)

# Plots
res <- r1
fits <- res$fits
reliab <- res$reliab
plt1_A <- plot_metrics(reliab, tols = res$confs_rel)
plt2_A <- plot_time_comparison_tol(fits, reliab, res$idx)
plt3_A <- plot_time_comparison_tol(fits, reliab, res$idx, TRUE)
diags_A <- get_diags_df(fits) # rhat and reff

# Plots
res <- r2
fits <- res$fits
reliab <- res$reliab
plt1_B <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_B <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_B <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_B <- get_diags_df(fits) # rhat and reff


# Plots
res <- r3
fits <- res$fits
reliab <- res$reliab
plt1_C <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_C <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_C <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_C <- get_diags_df(fits) # rhat and reff
