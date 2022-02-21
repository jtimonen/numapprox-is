#!/usr/bin/env Rscript
# TMDD experiment

# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Create model and simulation solver
solver_sim <- bdf(
  rel_tol = 1e-15,
  abs_tol = 1e-15,
  max_num_steps = 1e9
)
prior <- ode_model_tmdd(prior_only = TRUE)

# Define simulation parameters
sim_k <- c(0.592, 0.900, 2.212, 0.823, 0.201, 0.024)
sim_sigma <- 0.5
sim_params <- prior$make_params(c(sim_k, sim_sigma))

# Simulate ODE solution
t_sim <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, seq(2, 10, by = 1))
L0_sim <- 10
t0_sim <- 0
sim <- prior$gqs(
  t0 = t0_sim,
  t = t_sim,
  data = list(L0 = L0_sim, D = 3),
  params = sim_params,
  solver = solver_sim
)
y_sol <- sim$extract_odesol()
P_sol <- y_sol[1, , 3]

# Add noise and save data
SEED <- 123
set.seed(SEED)
P_dat <- P_sol + rnorm(P_sol, sd = sim_sigma)
dat <- list(
  t = t_sim, P_sol = P_sol, P_obs = P_dat, sim_k = sim_k,
  sim_sigma = sim_sigma, solver_sim = solver_sim, L0_sim = L0_sim,
  t0_sim = t0_sim, seed = SEED
)
saveRDS(dat, file = "simulated_data.rds")
