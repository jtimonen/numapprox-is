#!/usr/bin/env Rscript
library(rstan)

# Settings
theta    <- c(1.0, 2.0)
by       <- 0.25
T_max    <- 8
N_sets   <- 1
y0       <- c(1, 1)

# Compile simulation model
sm <- stan_model(file='stan/lv_simulate.stan')

# Seed
set.seed(123)
stan_seed <- 123
sigma <- 0.5

# Simulate data
th <- as.array(c(theta))
ts <- seq(by, T_max, by=by)
d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th, sigma=sigma)
f1 <- sampling(sm, data = d1, seed=stan_seed, chains = 1, cores = 1, iter = N_sets, warmup=0, algorithm="Fixed_param")
y_hat <- rstan::extract(f1)$y_hat
y  <- rstan::extract(f1)$y

# Save data
data_idx <- 1
dat <- list(y_hat=y_hat[data_idx,,], y=y[data_idx,,], 
            ts=ts, t0=0, T=length(ts))
fn <- paste0('data/dat_sigma_', sigma, '_set_', data_idx, '.rds')
saveRDS(dat, file=fn)
cat(paste0('Saved data to ', fn))

