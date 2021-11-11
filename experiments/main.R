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
  sig_figs = 18 # number of significant figures to store in floats
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
gpar <- paste("gamma[", c(1:4), "]", sep = "")
param_names <- c("beta", gpar, "phi_inv")
setup <- OdeExperimentSetup$new(
  "gsir", solver, solver_args_gen,
  stan_opts, param_names
)

# Fit prior model
prior_fit <- setup$sample_prior(refresh = 0)
prior_draws <- prior_fit$draws(setup$param_names)

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, prior_draws, setup$solver_args_gen)

# Plot generated data
setup$plot(prior_sim)
setup$add_simulated_data(prior_sim)
setup$set_init(prior_draws)
setup$plot(prior_sim)

# Run workflow
max_num_steps <- 1e4
run <- run_workflow(setup, 1e-4, 2, max_num_steps)

# Reference timing
tols <- 1 / run$tuning$metrics$inv_tol
tps <- setup$time_posterior_sampling(tols, max_num_steps, chains = 4)
t1_plot <- plot_timing(tols, tps$total)
t2_plot <- plot_timing(tols, tps$sampling)
t3_plot <- plot_timing(tols, tps$warmup)
