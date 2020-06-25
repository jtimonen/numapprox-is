#!/usr/bin/env Rscript
args     <- commandArgs(trailingOnly = TRUE)
DATA_IDX <- 1 #as.numeric(args[1])

# Setup
library(rstan)
library(loo)
source("functions.R")

# Compile model(s)
model_1 <- stan_model(file = 'stan/lv_rk45.stan')
#model_2 <- stan_model(file = 'stan/lv_rk4.stan')

# Options
ADAPT_DELTA <- 0.99
CHAINS      <- 4
ITER        <- 4000

# Load data
cat("-------------------------------------------------------------------------\n")
sigma <- 0.5
fn    <- paste0('data/dat_sigma_', sigma, '_set_', DATA_IDX, '.rds')
data  <- readRDS(file = fn)
cat(paste0('Read file ', fn, '\n'))

# Additional data
data$abs_tol_INF_  <- 1.0E-3
data$rel_tol_INF_  <- 1.0E-3
data$max_iter_INF_ <- 1.0E6
data$abs_tol_REF_  <- 1.0E-10
data$rel_tol_REF_  <- 1.0E-10
data$max_iter_REF_ <- 1.0E6

# Run inference 
new_data <- data
res <- run_inference(model_1, new_data, ITER, CHAINS, ADAPT_DELTA)
print(res$pareto_k)
print(res$runtimes)
