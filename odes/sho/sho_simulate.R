#!/usr/bin/env Rscript
library(rstan)

# Settings
col1     <- 'gray70'
col2     <- 'gray10'
lwd      <- 2
theta    <- 0.3
SIGMA    <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
by       <- 0.5
T_max    <- 10
N_sets   <- 30

# Compile simulation model
sm <- stan_model(file='sho_simulate.stan')

# Seed
set.seed(123)
stan_seed <- 123

for(i in 1:length(SIGMA)){
  sigma <- SIGMA[i]
  
  # Simulate data
  th <- as.array(c(theta))
  y0 <- c(1, -1)
  ts <- seq(by, T_max, by=by)
  d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th, sigma=sigma)
  f1 <- sampling(sm, data = d1, seed=stan_seed, chains = 1, cores = 1, iter = N_sets, warmup=0, algorithm="Fixed_param")
  y_hat <- rstan::extract(f1)$y_hat
  y  <- rstan::extract(f1)$y
  
  # Save data
  for(data_idx in 1:N_sets){
    dat <- list(y=y[data_idx,,], ts=ts, t0=0, T=length(ts))
    fn <- paste0('data/dat_sigma_', sigma, '_set_', data_idx, '.rds')
    saveRDS(dat, file=fn)
  }
}

# Plot data
#ylim <- c(min(y), max(y))
#xlim <- c(min(ts), max(ts))
#par(mfrow=c(2,1))
#plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[1]))
#lines(ts, y_hat[data_idx,,1], col=col1, lwd=lwd)
#points(ts, y[data_idx,,1], col=col2, pch=20)

#plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[2]))
#lines(ts, y_hat[data_idx,,2], col=col1, lwd=lwd)
#points(ts, y[data_idx,,2], col=col2, pch=20)

