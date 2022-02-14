#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results_rk45"
odemodeling:::create_dir_if_not_exist(res_dir)
source("setup.R")

# Run sampling
tols <- c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 10^(-5:-12))
solvers <- rk45_list(tols = tols, max_num_steps = 1e5)
fits <- model$sample_manyconf(
  t0 = dat$t0, t = dat$t, data = add_data, init = init,
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
