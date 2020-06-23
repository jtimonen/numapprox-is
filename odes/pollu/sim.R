#!/usr/bin/env Rscript

# Air pollution model from: https://archimede.dm.uniba.it/~testset/problems/pollu.php

# Settings
library(rstan)
col1     <- 'gray70'
col2     <- 'gray10'
lwd      <- 2
theta    <- c(0.35, 26.6, 12300, 0.00086, 0.00082, 
              15000, .00013, 24000, 16500, 9000,
              0.022, 12000, 1.88, 16300, 4800000,
              0.00035, 0.0175, 1e9, 0.444*1e12, 1240, 
              2.1, 5.78, 0.0474, 1780, 3.12)
SIGMA    <- c(0.01)
n_data   <- 60
T_max    <- 3000
N_sets   <- 1

# Compile simulation model
sm <- stan_model(file='stan/pollu_simulate.stan')

# Seed
set.seed(123)
stan_seed <- 123

for(i in 1:length(SIGMA)){
  sigma <- SIGMA[i]
  
  # Simulate data
  th <- as.array(c(theta))
  y0 <- c(1, -1)
  ts <- seq(by, T_max, length.out = n_data)
  d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th, sigma=sigma)
  f1 <- sampling(sm, data = d1, seed=stan_seed, chains = 1, cores = 1, iter = N_sets, warmup=0, algorithm="Fixed_param")
  y_hat <- rstan::extract(f1)$y_hat
  y  <- rstan::extract(f1)$y
  
  # Save data
  for(data_idx in 1:N_sets){
    dat <- list(y=y[data_idx,,], ts=ts, t0=0, T=length(ts), y_hat=y_hat[data_idx,,])
    fn <- paste0('data/dat_sigma_', sigma, '_set_', data_idx, '.rds')
    saveRDS(dat, file=fn)
  }
}
