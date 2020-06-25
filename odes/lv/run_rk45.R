#!/usr/bin/env Rscript
args     <- commandArgs(trailingOnly = TRUE)
DATA_IDX <- as.numeric(args[1])

# Setup
library(rstan)
library(loo)
source("functions.R")

# Compile model
model <- stan_model(file = 'stan/lv_rk45.stan')

# Options
ADAPT_DELTA <- 0.95
CHAINS      <- 100
ITER        <- 4000

# Load data
fn    <- paste0('data/dat_sigma_set_', DATA_IDX, '.rds')
data  <- readRDS(file = fn)
cat(paste0('Read data from file ', fn, '\n'))
data$abs_tol_REF_  <- 1.0E-10
data$rel_tol_REF_  <- 1.0E-10
data$max_iter_REF_ <- 1.0E6

# Run with different tolerances
TOL <- 10^(-10:-3)
TOL <- c(TOL, 0.005)

L <- length(TOL)
for(i in 1:L){
  tol <- TOL[i]
  
  # Additional data
  new_data <- data
  new_data$abs_tol_INF_  <- tol
  new_data$rel_tol_INF_  <- tol
  new_data$max_iter_INF_ <- 1.0E6
  
  # Run inference 
  res <- run_inference(model, new_data, ITER, CHAINS, ADAPT_DELTA)
  print(res$pareto_k)
  print(res$runtimes)
  res$tol <- tol
  fn <- paste0('res/res_rk45_', i ,'.rds')
  cat(paste0('Saving to file  ', fn, '\n'))
  saveRDS(res, fn)
}

