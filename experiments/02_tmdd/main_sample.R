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
source("setup_tmdd.R")

# Create experiment setup
solver_args_gen <- list(
  rel_tol = 1e-9,
  abs_tol = 1e-9,
  max_num_steps = 1e9
)
solver <- "rk45"
kpar <- paste("k[", c(1:6), "]", sep = "")
param_names <- c(kpar, "sigma")
setup <- OdeExperimentSetup$new(
  "tmdd", solver, solver_args_gen,
  stan_opts, param_names
)
setup$set_hmc_initial_step_size(0.1)
print(setup)

# Fit prior model
prior_fit <- setup$sample_prior(refresh = 0)
prior_draws <- prior_fit$draws(setup$param_names)

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, prior_draws, setup$solver_args_gen)

# Plot and save generated data
setup$add_simulated_data(prior_sim)
setup$set_init(prior_draws)
plot_prior <- setup$plot(prior_sim)
cat("SETUP:\n")
print(setup)
cat("DATA:\n")
print(setup$data)
cat("INIT:\n")
print(setup$init)
saveRDS(setup$data, file = fn_data)

# Run workflow
if (idx > 40) {
  MNS <- 1e5
} else if (idx > 20) {
  MNS <- 1e4
} else {
  MNS <- 1e3
}
max_num_steps <- MNS
run <- run_workflow(setup, 1e-3, 10, max_num_steps, 7)

# Save result
run$post_fit$save_object(fn_fit)

# Reference timing
# tols <- 1 / run$tuning$metrics$inv_tol
tols <- 10^(c(-4, -6, -8, -10, -12))
tps <- setup$time_posterior_sampling(res_dir, idx, tols, max_num_steps, chains = 4)

seed <- run$post_fit$runset$args$seed
all_results <- list(
  run = run, tps = tps, plot_prior = plot_prior,
  max_num_steps = max_num_steps, setup = setup,
  seed = seed
)
saveRDS(all_results, file = fn_res)
siz <- format(object.size(all_results), units = "Kb")
cat("save size:", siz, "\n")
