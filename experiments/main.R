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

# Options
stan_opts <- list(
  sig_figs = 18 # number of significant figures to store in floats
)

# R functions
source("functions.R")
source("setup_sir.R")

# Create experiment setup
solver_args_gen <- list(
  rel_tol = 1e-9,
  abs_tol = 1e-9,
  max_num_steps = 1e9
)
solver <- "rk45"
param_names <- c("beta", "gamma", "phi_inv")
setup <- OdeExperimentSetup$new(
  "sir", solver, solver_args_gen,
  stan_opts, param_names
)

# Fit prior model
prior_fit <- setup$sample_prior(refresh = 1000)
prior_draws <- prior_fit$draws(setup$param_names)

# Simulate solutions and data using prior draws
prior_sim <- simulate(setup, prior_draws, setup$solver_args_gen)

# Plot generated data
setup$plot(prior_sim)
setup$add_simulated_data(prior_sim)
setup$set_init(prior_draws)
setup$plot(prior_sim)

# Sample posterior
TOLS <- 10^seq(-4,-8)
L <- length(TOLS)
WT <- matrix(0.0, L, 4)
ST <- matrix(0.0, L, 4)
TT <- matrix(0.0, L, 4)
j <- 0
for (ttol in TOLS) {
  print(ttol)
  j <- j + 1
  solver_args_sample <- list(
    abs_tol = ttol, #1e-4,
    rel_tol = ttol, #1e-4,
    max_num_steps = 1e6
  )
  post_fit <- setup$sample_posterior(solver_args_sample, refresh = 1000)
  WT[j,] <- post_fit$time()$chains$warmup
  ST[j,] <- post_fit$time()$chains$sampling
  TT[j,] <- post_fit$time()$chains$total
}

post_draws <- post_fit$draws(setup$param_names)

# Simulate using posterior draws
post_sim <- simulate(setup, post_draws, solver_args_sample)
plt_sim <- setup$plot(post_sim)

# Tune the numerical method in $M_{high}$ so that it is reliable at these draws.
tuning_factor <- 2
tuning <- tune_solver_tols(
  setup, post_sim, post_draws, solver_args_sample, tuning_factor
)
plt_tuning <- plot_tuning(tuning, nrow=1)

plt <- ggarrange(plt_tuning, plt_sim, nrow=2, labels=c("", "e"))

# Compute importance weights and Pareto-$k$
is <- use_psis(post_sim, tuning$last_sim)

# Importance resampling
post_draws_resampled <- posterior::resample_draws(post_draws, weights = is$log_weights)

# Comparing runtimes
workflow_time <- post_fit$time()$total + post_sim$time()$total + tuning$total_time
cat("Total workflow time was", workflow_time, "seconds.\n", sep = " ")


solver_args_refsample <- list(
  abs_tol = tuning$last_tol,
  rel_tol = tuning$last_tol,
  max_num_steps = solver_args_sample$max_num_steps
)
fit_ref <- sample_posterior(stanmodels$posterior, dat,
                            solver_args_refsample, stan_opts,
                            refresh = 1000,
                            init = 0
)
print(fit_ref$time()$total)
