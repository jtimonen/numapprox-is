#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results_bdf"
odemodeling:::create_dir_if_not_exist(res_dir)
source("setup.R")

# Run sampling
tols <- seq(0.06, 0.01, by = -0.01)
tols <- c(tols, 0.001, 0.0001, 10^(-5:-12))
solvers <- bdf_list(tols = tols, max_num_steps = 1e5)
fits <- model$sample_manyconf(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solvers = solvers, step_size = 0.1, savedir = res_dir,
  iter_warmup = ITER, iter_sampling = ITER, chains = CHAINS
)

# Save results
results <- list(
  fits = fits,
  session_info = sessionInfo()
)
fp <- file.path(res_dir, "sampling.rds")
saveRDS(results, file = fp)
