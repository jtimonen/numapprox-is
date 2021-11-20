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
fn <- file.path(res_dir, paste0("res_", idx, ".rds"))
cat("Results will be saved to: ", fn, "\n", sep = "")
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
source("classes.R")
source("functions.R")
source("setup_gsir.R")

# Create experiment setup
solver_args_gen <- list(
  rel_tol = 1e-9,
  abs_tol = 1e-9,
  max_num_steps = 1e9
)
solver <- "rk45"
gpar <- paste("gamma[", c(1:10), "]", sep = "")
param_names <- c("beta", gpar, "phi_inv")
setup <- OdeExperimentSetup$new(
  "gsir", solver, solver_args_gen,
  stan_opts, param_names
)
print(setup)

# Fit prior model
prior_fit <- setup$sample_prior(refresh = 0)
prior_draws <- prior_fit$draws(setup$param_names)

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, prior_draws, setup$solver_args_gen)

# Plot generated data
setup$add_simulated_data(prior_sim)
setup$set_init(prior_draws)
plot_prior <- setup$plot(prior_sim)
cat("SETUP:\n")
print(setup)
cat("DATA:\n")
print(setup$data)
cat("INIT:\n")
print(setup$init)

# Run workflow
if (idx > 40) {
  MNS <- 1e6
} else if (idx > 20) {
  MNS <- 1e5
} else {
  MNS <- 1e4
}
MNS <- 1e6
max_num_steps <- MNS
run <- run_workflow(setup, 1e-3, 10, max_num_steps, 6)

# Reference timing
# tols <- 1 / run$tuning$metrics$inv_tol
tols <- c(1e-4, 1e-6)
tps <- setup$time_posterior_sampling(tols, max_num_steps, chains = 4)
t1_plot <- plot_timing(tols, tps$total)
t2_plot <- plot_timing(tols, tps$sampling)
t3_plot <- plot_timing(tols, tps$warmup)

seed <- run$post_fit$runset$args$seed
all_results <- list(
  run = run, tps = tps, plot_prior = plot_prior,
  max_num_steps = max_num_steps, setup = setup,
  seed = seed
)
siz <- format(object.size(all_results), units = "Kb")
cat("save size:", siz, "\n")
saveRDS(all_results, file = fn)
