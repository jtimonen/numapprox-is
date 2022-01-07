#!/usr/bin/env Rscript
# args     <- commandArgs(trailingOnly = TRUE)
# DATA_IDX <- as.numeric(args[1])

# Setup
library(rstan)
library(loo)
source("functions.R")

# Compile model
model <- stan_model(file = "stan/lv_rk4.stan")

# Options
ADAPT_DELTA <- 0.95
CHAINS <- 30
ITER <- 4000
for (DATA_IDX in 1:10) {

  # Load data
  fn <- paste0("data/dat_set_", DATA_IDX, ".rds")
  data <- readRDS(file = fn)
  cat(paste0("Read file ", fn, "\n"))
  data$abs_tol_REF_ <- 1.0E-10
  data$rel_tol_REF_ <- 1.0E-10
  data$max_iter_REF_ <- 1.0E6

  # Run with different step sizes
  H <- seq(0.1, 1, by = 0.1)
  L <- length(H)
  for (i in 1:L) {
    step_size <- H[i]

    # Additional data
    new_data <- data
    new_data <- add_interpolation_data(new_data, h = step_size)

    # Run inference
    res <- run_inference(model, new_data, ITER, CHAINS, ADAPT_DELTA)
    print(res$pareto_k)
    print(res$runtimes)
    res$step_size <- step_size
    fn <- paste0("res/rk4/rk4_dat_", DATA_IDX, "_step_", i, ".rds")
    cat(paste0("Saving to file  ", fn, "\n"))
    saveRDS(res, fn)
  }
}
