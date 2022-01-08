#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# R functions and requirements
source("../utils.R")
source("../models.R")
library(odemodeling)
library(posterior)

# Setup
idx <- setup_experiment_index(args)
fp <- setup_experiment_paths(idx)

# Options
stan_opts <- list(
  sig_figs = 12, # number of significant figures to store in floats
  seed = 123
)

# Create experiment setup
solver_gen <- bdf(
  rel_tol = 1e-15,
  abs_tol = 1e-15,
  max_num_steps = 1e9
)
prior <- ode_model_tmdd(prior_only = TRUE)

# SIMULATION ----------------------------------------------------------

# Define simulation parameters
param_names <- c("k_on", "k_off", "k_in", "k_out", "k_eP", "k_eL", "sigma")
sim_k <- c(0.592, 0.900, 2.212, 0.823, 0.201, 0.024)
sim_sigma <- 0.3
sim_params <- as_draws_array(array(c(sim_k, sim_sigma), dim = c(1, 1, 7)))
dimnames(sim_params)$variable <- param_names

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, sim_params, setup$solver_args_gen)

# Plot and save generated data
setup$add_simulated_data(prior_sim)
setup$set_init(init = 0)
plot_prior <- setup$plot(prior_sim)
print(setup)
saveRDS(setup$data, file = fn_data)

# Denser plot
h <- 0.01
new_t <- seq(h, max(setup$data$t), by = h)
sim_dense <- simulate_dense(setup, sim_params, setup$solver_args_gen, new_t)
x_dense <- get_x_sim(sim_dense)
df_dense <- data.frame(x_dense, rep(c("L", "R", "P"), each = 2000), rep(new_t, times = 3))
colnames(df_dense) <- c("x", "var", "t")
plt_sim <- ggplot(df_dense, aes(x = t, y = x, group = var, color = var)) +
  geom_line(lwd = 1) +
  xlab("Concentration") +
  ylab("Time")
df_data <- data.frame(setup$data$t, setup$data$y)
colnames(df_data) <- c("t", "y")
plt_sim <- plt_sim + geom_point(
  data = df_data, aes(x = t, y = y),
  inherit.aes = FALSE
) + ggtitle("Black dots = data of L (drug concentration)")
ggsave("data_tmdd.pdf", width = 8, height = 5.3)

# SAMPLING ----------------------------------------------------------
max_num_steps <- 1e6
tols <- c(
  0.1, 0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
  1e-11, 1e-12
)
post <- setup$sample_posterior_many(res_dir, idx, tols, max_num_steps, chains = 4)

# Run workflow
idx <- 2
sampled <- load_fit(post, idx)
L <- length(tols)
tols_val <- tols[(idx + 1):L]
run <- validate_fit(setup, sampled, tols_val, max_num_steps)

# Save result
tune <- run$tuning
tp <- plot_tuning(tune)
ggsave("tmdd_tuning.pdf", width = 9.67, height = 6.17)

# Plot times
t <- post$grand_total[idx:L]
tols <- post$tols[idx:L]
plot(-log10(tols), t, "o",
  ylab = "time (s)", pch = 16, ylim = c(0, 320),
  xlab = "T", xaxt = "n"
)
grid()
t_sample <- t[1]
t2 <- tune$time + t_sample
tols2 <- 1 / tune$inv_tol
lines(-log10(tols2), t2, col = "firebrick3")
points(-log10(tols2), t2, col = "firebrick3", pch = 17)
legend(2, 200, c("HMC-NUTS using tol=T", "HMC-NUTS using tol=0.05 + PSIS with tol=T"),
  lty = c(1, 1), col = c("black", "firebrick3"),
  pch = c(16, 17)
)
axis(1, at = -log10(tols), las = 2, labels = tols)
