#!/usr/bin/env Rscript

# Settings
library(rstan)
col1     <- 'gray70'
col2     <- 'gray10'
lwd      <- 2

# Read data
dat_real <- readRDS('data/data13.rds')
dat_sim <- list()

# Setup
dat_sim$K        <- 9
dat_sim$G        <- 60
dat_sim$D        <- 42
dat_sim$S        <- 42
dat_sim$ts       <- c(1:42)
dat_sim$tswitch  <- 20
dat_sim$age_dist <- dat_real$age_dist
dat_sim$pop_t    <- dat_real$pop_t
dat_sim$contact  <- dat_real$contact
dat_sim$p_gamma  <- dat_real$p_gamma

SIGMA    <- c(0.1)
by       <- 0.05
T_max    <- 3
N_sets   <- 1

# BDF options
rtol <- 1e-10
atol <- 1e-10
max_steps <- 1e6

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
  ts <- seq(by, T_max, by = by)
  d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th, sigma=sigma,
             tol1=rtol, tol2=atol, max_steps=max_steps)
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

n_data <- length(ts)
y_last <- y_hat[1,n_data,]
print(y_last - y_ref)

#y_ref <- y_last


