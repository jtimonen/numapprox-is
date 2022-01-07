#!/usr/bin/env Rscript
library(rstan)

# Settings
theta <- c(1.0, 2.0)
by <- 0.5
T_max <- 6
N_sets <- 10
y0 <- c(1, 1)

# Compile simulation model
sm <- stan_model(file = "stan/lv_simulate.stan")

# Seed
sigma <- 0.5

# Simulate data
th <- as.array(c(theta))
ts <- seq(by, T_max, by = by)
d1 <- list(T = length(ts), y0 = y0, t0 = 0, ts = ts, theta = th, sigma = sigma)
f1 <- sampling(sm, data = d1, chains = 1, cores = 1, iter = N_sets, warmup = 0, algorithm = "Fixed_param")
y_hat <- rstan::extract(f1)$y_hat
y <- rstan::extract(f1)$y

# Save data
for (data_idx in 1:N_sets) {
  dat <- list(
    y_hat = y_hat[data_idx, , ], y = y[data_idx, , ],
    ts = ts, t0 = 0, T = length(ts), sigma = sigma, idx = data_idx
  )
  fn <- paste0("data/dat_set_", data_idx, ".rds")
  saveRDS(dat, file = fn)
  cat(paste0("Saved data to ", fn, "\n"))
}
