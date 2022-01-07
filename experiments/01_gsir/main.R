#!/usr/bin/env Rscript
if (interactive()) {
  idx <- 0
} else {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- args[1]
}
res_dir <- "res"
if (!dir.exists(res_dir)) {
  message("res_dir doesn't exist, creating it...")
  dir.create(res_dir)
}

cat("\n------ idx = ", idx, " --------\n", sep = "")
fn_res <- file.path(res_dir, paste0("res_", idx, ".rds"))
fn_fit <- file.path(res_dir, paste0("fit_", idx, ".rds"))
fn_data <- file.path(res_dir, paste0("dat_", idx, ".rds"))
cat("Results will be saved to: ", fn_res, "\n", sep = "")
cat("Fit will be saved to: ", fn_fit, "\n", sep = "")
cat("Data will be saved to: ", fn_data, "\n", sep = "")
idx <- as.numeric(idx)

# Requirements
library(cmdstanr)
library(posterior)
library(bayesplot)
library(checkmate)
library(loo)
library(stats)
library(outbreaks)
library(scales)
library(ggplot2)
library(ggdist)
library(R6)
library(ggpubr)
library(tidyr)

# Options
stan_opts <- list(
  sig_figs = 12 # number of significant figures to store in floats
)

# R functions
source("../R/classes.R")
source("../R/functions.R")
source("setup_gsir.R")

# Create experiment setup
solver_args_gen <- list(
  rel_tol = 1e-13,
  abs_tol = 1e-13,
  max_num_steps = 1e9
)
solver <- "rk45"
G <- 8
gpar <- paste("gamma[", c(1:G), "]", sep = "")
ppar <- paste("phi_inv[", c(1:G), "]", sep = "")
param_names <- c("beta", gpar, ppar)
setup <- OdeExperimentSetup$new(
  "gsir", solver, solver_args_gen,
  stan_opts, param_names
)
setup$set_hmc_initial_step_size(0.1)
print(setup)

# Run workflow
if (idx > 40) {
  MNS <- 1e5
} else if (idx > 20) {
  MNS <- 1e4
} else {
  MNS <- 1e3
}


# SIMULATION ----------------------------------------------------------

# Fit prior model
prior_fit <- setup$sample_prior(refresh = 0)
prior_draws <- prior_fit$draws(setup$param_names)

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, prior_draws, setup$solver_args_gen)

# Plot and save generated data
setup$add_simulated_data(prior_sim)
setup$set_init(init = 0)
plot_prior <- setup$plot(prior_sim)
print(setup)
saveRDS(setup$data, file = fn_data)

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
