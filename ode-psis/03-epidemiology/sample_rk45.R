#!/usr/bin/env Rscript

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Setup data and model
source("setup.R")
ITER <- 1000
res_dir <- "results_rk45"
odemodeling:::create_dir_if_not_exist(res_dir)

# Sampling
tols <- c(0.05, 0.001, 10^c(-2:-10))
solvers <- rk45_list(tols = tols, max_num_steps = 1e3)
fits <- model$sample_manyconf(
  t0 = t0, t = t, data = add_data, init = init,
  solvers = solvers,
  step_size = step_size, iter_warmup = ITER, iter_sampling = ITER, chains = 4,
  adapt_delta = 0.9, max_treedepth = 13, refresh = 25,
  savedir = res_dir
)

# Save results
results <- list(
  fits = fits,
  session_info = sessionInfo()
)
fp <- file.path(res_dir, "sampling.rds")
saveRDS(results, file = fp)
