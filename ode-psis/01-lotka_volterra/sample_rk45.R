#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results_rk45"
odemodeling:::create_dir_if_not_exist(res_dir)

# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Create model and load data
model <- ode_model_lv()
dat <- load_data_lynxhare()
add_data <- list(
  y_obs_init = dat$y_obs_init,
  y_obs = dat$y_obs,
  D = 2
)

# Define initial params for sampling
init <- list(
  alpha = 1,
  beta = 0.1,
  gamma = 1,
  delta = 0.1,
  y0 = dat$y_obs_init,
  sigma = c(1, 1)
)
init <- rep(list(init), CHAINS) # same for all chains

# Run sampling
tols <- c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 10^(-5:-12))
solvers <- rk45_list(tols = tols, max_num_steps = 1e6)
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
