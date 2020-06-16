#!/usr/bin/env Rscript
args     <- commandArgs(trailingOnly = TRUE)
DATA_IDX <- as.numeric(args[1])

# Setup
library(rstan)
library(loo)
source("functions.R")
model_1 <- stan_model(file = 'stan/sho_rk45.stan')
model_2 <- stan_model(file = 'stan/sho_rk4.stan')

# Options
ADAPT_DELTA <- 0.8
CHAINS      <- 4
ITER        <- 2000
SIGMA       <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
STEP_SIZE   <- c(0.05, 0.1, 0.2, 0.5, 1.0)
S           <- length(SIGMA)
J           <- length(STEP_SIZE)
PARETO_K    <- matrix(0, S, J+1)
RUNTIMES    <- array(0, c(S, J+1, CHAINS))

for(i in 1:S){
  
  # Load data
  cat("-------------------------------------------------------------------------\n")
  sigma <- SIGMA[i]
  fn    <- paste0('data/dat_sigma_', sigma, '_set_', DATA_IDX, '.rds')
  data  <- readRDS(file = fn)
  cat(paste0('Read file ', fn, '\n'))
  
# Additional data for reference method ode_integrate_rk45
  data$abs_tol_REF_  <- 1.0E-6
  data$rel_tol_REF_  <- 1.0E-6
  data$max_iter_REF_ <- 1.0E3
  
  # Run inference with the reference model
  res <- run_inference(model_1, data, ITER, CHAINS, ADAPT_DELTA)
  print(res)

  # Run inference with different step sizes of the approximate model
  for(j in 1:J){
      new_data  <- data
      step_size <- STEP_SIZE[j]
      new_data  <- add_interpolation_data(new_data, step_size)
      res       <- run_inference(model_2, new_data, ITER, CHAINS, ADAPT_DELTA)
      print(j)
      print(res)
  }

}

run_inference <- function(model, data, ITER, CHAINS, ADAPT_DELTA){

  # Run sampling
  fit <- sampling(object  = model,
                  data    = data,
                  iter    = ITER,
                  chains  = CHAINS,
                  control = list(adapt_delta = ADAPT_DELTA))
  
  # Extract log posterior values (not Jacobian adjusted)
  lh1 <- get_samples(fit, 'log_lik_na')
  lh2 <- get_samples(fit, 'log_lik_na_REF_')
  pr1 <- get_samples(fit, 'log_prior_na')
  pr2 <- get_samples(fit, 'log_prior_na_REF_')
  post1 <- lh1 + pr1
  post2 <- lh2 + pr2
  
  # PSIS
  out <- psis(post1 - post2)
  pareto_k <- out$diagnostics$pareto_k
  cat(paste0('>>> sigma=', sigma, ', pareto_k=', pareto_k, '\n'))
  runtimes <- rowSums(get_elapsed_time(fit))
  
  # Return list
  return(list(pareto_k=pareto_k, runtimes=runtimes))

}

print(RUNTIMES)
print(PARETO_K)

