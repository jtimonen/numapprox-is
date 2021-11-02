#!/usr/bin/env Rscript

# Settings
library(rstan)
by <- 1
T_max <- 40
N_sets <- 1
gamma <- 0.2
beta <- 0.7
phi_inv <- 0.2

# Solver options
rtol <- 1e-10
atol <- 1e-10
max_iter <- 1e6

# Compile simulation model
sm <- stan_model(file = "stan/sir_simulate.stan")

# Seed
set.seed(12)
stan_seed <- 12

# Simulate data
data_idx <- 1
ts <- seq(by, T_max, by = by)
opt <- list(
  T = length(ts), t0 = 0, ts = ts,
  gamma = gamma, beta = beta, phi_inv = phi_inv,
  abs_tol = atol, rel_tol = rtol, max_iter = max_iter
)
sim <- sampling(sm,
  data = opt, seed = stan_seed, chains = 1, cores = 1,
  iter = N_sets, warmup = 0, algorithm = "Fixed_param"
)
y_hat <- rstan::extract(sim)$y_hat
y <- rstan::extract(sim)$y

# Save data
for (data_idx in 1:N_sets) {
  yi <- y[data_idx, ]
  dat <- list(
    ts = ts, t0 = 0, T = length(ts),
    y = t(t(yi)), y_hat = y_hat[data_idx, , ],
    beta = beta, gamma = gamma, phi_inv = phi_inv
  )
  fn <- paste0("data/dat_set_", data_idx, ".rds")
  saveRDS(dat, file = fn)
}

n_data <- length(ts)
y_last <- y_hat[1, n_data, ]
# print(y_last - y_ref)

# y_ref <- y_last
