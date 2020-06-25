#!/usr/bin/env Rscript
args     <- commandArgs(trailingOnly = TRUE)
DATA_IDX <- 1 #as.numeric(args[1])

# Setup
library(rstan)
library(loo)
source("functions.R")

# Compile model(s)
model_1 <- stan_model(file = 'stan/lv_rk45.stan')

# Options
ADAPT_DELTA <- 0.95
CHAINS      <- 4
ITER        <- 2000

# Load data
cat("-------------------------------------------------------------------------\n")
sigma <- 0.5
fn    <- paste0('data/dat_sigma_', sigma, '_set_', DATA_IDX, '.rds')
data  <- readRDS(file = fn)
cat(paste0('Read file ', fn, '\n'))
data$abs_tol_REF_  <- 1.0E-10
data$rel_tol_REF_  <- 1.0E-10
data$max_iter_REF_ <- 1.0E6

# Additional data
new_data <- data
new_data$abs_tol_INF_  <- 1.0E-6
new_data$rel_tol_INF_  <- 1.0E-6
new_data$max_iter_INF_ <- 1.0E6

# Run inference 
new_data <- data
res_1 <- run_inference(model_1, new_data, ITER, CHAINS, ADAPT_DELTA)
print(res_1$pareto_k)
print(res_1$runtimes)

# Compile other model
model_2 <- stan_model(file = 'stan/lv_rk4.stan')

# Additional data
new_data <- data
new_data <- add_interpolation_data(new_data, h = 0.2)

# Run inference 
res_2 <- run_inference(model_2, new_data, ITER, CHAINS, ADAPT_DELTA)
print(res_2$pareto_k)
print(res_2$runtimes)
