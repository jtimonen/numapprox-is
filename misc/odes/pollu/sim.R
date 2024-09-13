#!/usr/bin/env Rscript

# Air pollution model from: https://archimede.dm.uniba.it/~testset/problems/pollu.php

# Settings
library(rstan)
col1 <- "gray70"
col2 <- "gray10"
lwd <- 2
k <- rep(0, 25)

k[1] <- 0.350 * 1e0
k[2] <- 0.266 * 1e2
k[3] <- 0.123 * 1e5
k[4] <- 0.860 * 1e-3
k[5] <- 0.820 * 1e-3
k[6] <- 0.150 * 1e5
k[7] <- 0.130 * 1e-3
k[8] <- 0.240 * 1e5
k[9] <- 0.165 * 1e5

k[10] <- 0.900 * 1e4
k[11] <- 0.220 * 1e-1
k[12] <- 0.120 * 1e5
k[13] <- 0.188 * 1e1
k[14] <- 0.163 * 1e5
k[15] <- 0.480 * 1e7
k[16] <- 0.350 * 1e-3
k[17] <- 0.175 * 1e-1
k[18] <- 0.100 * 1e9

k[19] <- 0.444 * 1e12
k[20] <- 0.124 * 1e4
k[21] <- 0.210 * 1e1
k[22] <- 0.578 * 1e1
k[23] <- 0.474 * 1e-1
k[24] <- 0.178 * 1e4
k[25] <- 0.312 * 1e1


SIGMA <- c(0.1)
by <- 0.1
T_max <- 3
N_sets <- 1

# BDF options
rtol <- 1e-10
atol <- 1e-10
max_steps <- 1e6

# Compile simulation model
sm <- stan_model(file = "stan/pollu_simulate.stan")

# Seed
set.seed(123)
stan_seed <- 123

for (i in 1:length(SIGMA)) {
  # Simulate data
  sigma <- SIGMA[i]
  theta <- as.array(k)
  y0 <- c(1, -1)
  ts <- seq(by, T_max, by = by)
  d1 <- list(
    T = length(ts), y0 = y0, t0 = 0, ts = ts, theta = theta, sigma = sigma,
    tol1 = rtol, tol2 = atol, max_steps = max_steps
  )
  f1 <- sampling(sm, data = d1, seed = stan_seed, chains = 1, cores = 1, iter = N_sets, warmup = 0, algorithm = "Fixed_param")
  y_hat <- rstan::extract(f1)$y_hat
  y <- rstan::extract(f1)$y

  # Save data
  for (data_idx in 1:N_sets) {
    dat <- list(
      y = y[data_idx, , ], ts = ts, t0 = 0, T = length(ts),
      y_hat = y_hat[data_idx, , ], k_data = k, sigma_data = sigma
    )
    fn <- paste0("data/dat_sigma_", sigma, "_set_", data_idx, ".rds")
    saveRDS(dat, file = fn)
  }
}

n_data <- length(ts)
y_last <- y_hat[1, n_data, ]
# print(y_last - y_ref)

# y_ref <- y_last
