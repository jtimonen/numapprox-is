#!/usr/bin/env Rscript
args     <- commandArgs(trailingOnly = TRUE)
#DATA_IDX <- as.numeric(args[1])
DATA_IDX <- 1

# Setup
library(rstan)
library(loo)
source("functions.R")

# Compile models
model_1 <- stan_model(file = 'stan/sir_rk45.stan')
#model_2 <- stan_model(file = 'stan/sho_rk4.stan')

# Options
ADAPT_DELTA <- 0.9
CHAINS      <- 4
ITER        <- 2000

# Load data
cat("-------------------------------------------------------------------------\n")
fn    <- paste0('data/dat_set_', DATA_IDX, '.rds')
data  <- readRDS(file = fn)
cat(paste0('Read file ', fn, '\n'))

# Additional data
data$abs_tol_REF_  <- 1.0E-10
data$rel_tol_REF_  <- 1.0E-10
data$max_iter_REF_ <- 1.0E6

# Run inference
new_data <- data
new_data$abs_tol_INF_  <- 1.0E-2
new_data$rel_tol_INF_  <- 1.0E-2
new_data$max_iter_INF_ <- 1.0E3
res <- run_inference(model_1, new_data, ITER, CHAINS, ADAPT_DELTA)

# Save solutions
YHAT <- get_samples(res$fit, 'y_hat')
YHAT_REF_ <- get_samples(res$fit, 'y_hat_REF_')
saveRDS(list(inf=YHAT, ref=YHAT_REF_), 'data/yhat.rds')

