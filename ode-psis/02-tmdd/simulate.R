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
solver_gen <- bdf(
  rel_tol = 1e-15,
  abs_tol = 1e-15,
  max_num_steps = 1e9
)
prior <- ode_model_tmdd(prior_only = TRUE)

# Define simulation parameters
sim_k <- c(0.592, 0.900, 2.212, 0.823, 0.201, 0.024)
sim_sigma <- 0.3
sim_params <- prior$make_params(c(sim_k, sim_sigma))

# Simulate ODE solution
t_sim <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, seq(2, 10, by = 1))
sim <- prior$gqs(
  t0 = 0,
  t = t_sim,
  data = list(L0 = 10, D = 3),
  params = sim_params,
  solver = solver_gen
)
y_sol <- sim$extract_odesol()
P_sol <- y_sol[1, , 3]

# Add noise and save data
P_dat <- P_sol + rnorm(P_sol, sd = sim_sigma)
dat <- list(
  t = t_sim, P_sol = P_sol, P_obs = P_dat, sim_k = sim_k,
  sim_sigma = sim_sigma
)
saveRDS(dat, file = "simulated_data.rds")

# Plot
plt_sim <- sim$plot_odesol(
  ydim_names = c("L", "R", "P"),
  include_y0 = TRUE
)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep("P", length(t_sim)))
plt <- plt_sim + 
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab("Concentration") + xlab("Time")

# Save plot
odemodeling:::create_dir_if_not_exist("figures")
fp  <- file.path("figures", "simulated_data.pdf")
ggsave(plt, file = fp, width = 7, height = 3.5)
