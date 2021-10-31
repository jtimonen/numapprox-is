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

# Options
stan_opts <- list(
  seed = 2353, # random seed for Stan
  sig_figs = 18 # number of significant figures to store in floats
)

# R functions
source("functions.R")
source("setup_sir.R")

data_list <- setup_standata_sir()
code_parts <- setup_stancode_sir(solver = "rk45")
stanmodels <- create_cmdstan_models(code_parts)

prior_fit <- stanmodels$prior$sample(
  sig_figs = stan_opts$sig_figs,
  seed = stan_opts$seed,
  refresh = 1000
)
print(prior_fit)
prior_draws <- prior_fit$draws(c("beta", "gamma", "phi_inv"))

solver_args_default <- list(
  rel_tol = 1e-6,
  abs_tol = 1e-6,
  max_num_steps = 1e6
)
# Simulate solutions using prior draws
prior_sim <- simulate(
  stanmodels$sim, prior_draws, data_list,
  solver_args_default, stan_opts
)

# Plot solutions

# Function that plots SIR solutions against data (of infected)
plot_sir_example_solutions <- function(sim, data, thin = 1, main = "") {
  N <- data$N
  x_sim <- posterior::thin_draws(posterior::merge_chains(sim$draws("x")), thin)
  I_sim <- x_sim[, 1, (N + 1):(2 * N), drop = TRUE]
  plot(data$t, rep(data$pop_size, data$N),
    type = "l", lty = 2,
    col = "gray70", ylim = c(0, 800), ylab = "Infected", xlab = "Day",
    main = main
  )
  for (i_draw in 1:nrow(I_sim)) {
    I <- as.vector(I_sim[i_draw, ])
    lines(data$t, I, col = scales::alpha("firebrick", 0.1))
  }
  points(data$t, data$I_data, ylim = c(0, 1000), pch = 20)
}

title <- "Solutions using prior param draws"
plot_sir_example_solutions(prior_sim, data_list, thin = 10, main = title)

## Analyzing the numerical method
TOLS <- 10^seq(-9, -4)
atols <- TOLS
rtols <- TOLS
sims <- simulate_many(
  stanmodels$sim, prior_draws, data_list, stan_opts,
  atols, rtols, solver_args_default$max_num_steps
)

# Plot
mean_abs_sol_error <- compute_sol_errors(sims$sims, "mean")
max_abs_sol_error <- compute_sol_errors(sims$sims, "max")
plot_sim_errors(atols, rtols, mean_abs_sol_error)
plot_sim_errors(atols, rtols, max_abs_sol_error)

mean_abs_loglik_error <- compute_loglik_errors(sims$log_liks, "mean")
max_abs_loglik_error <- compute_loglik_errors(sims$log_liks, "max")
plot_sim_errors(atols, rtols, mean_abs_loglik_error, log = FALSE)
plot_sim_errors(atols, rtols, max_abs_loglik_error, log = FALSE)

## Computational challenges
plot_sim_times(atols, rtols, sims$times)

solver_args_sample <- list(
  abs_tol = 1e-4,
  rel_tol = 1e-4,
  max_num_steps = 1e6
)
fit <- sample_posterior(stanmodels$posterior, data_list,
  solver_args_sample, stan_opts,
  refresh = 1000,
  init = 0
)

print(fit$time())
print(fit)
post_draws <- fit$draws(c("beta", "gamma", "phi_inv"))

# Simulate and plot generated solutions
post_sim <- simulate(stanmodels$sim, post_draws, data_list, solver_args_sample, stan_opts)

title <- "Solutions using posterior param draws"
plot_sir_example_solutions(post_sim, data_list, thin = 10, main = title)

# Tune the numerical method in $M_{high}$ so that it is reliable at these draws.
TOLS <- 10^seq(-4, -12) # could be just halving for example
tuning_tol <- 0.0001
tuning <- tune_solver(
  TOLS, stanmodels$sim, post_draws, data_list, stan_opts,
  solver_args_sample$max_num_steps, tuning_tol
)

# Print the final tolerance
print(tuning$last_tol)

# Compute importance weights and Pareto-$k$
is <- use_psis(post_sim, tuning$last_sim)
print(is)
print(is$diagnostics$pareto_k)

# Importance resampling
post_draws_resampled <- posterior::resample_draws(post_draws, weights = is$log_weights)

# Comparing runtimes
workflow_time <- fit$time()$total + sum(tuning$times)
cat("Total workflow time was", workflow_time, "seconds.\n", sep = " ")
solver_args_refsample <- list(
  abs_tol = tuning$last_tol,
  rel_tol = tuning$last_tol,
  max_num_steps = solver_args_sample$max_num_steps
)
fit_ref <- sample_posterior(stanmodels$posterior, data_list,
  solver_args_refsample, stan_opts,
  refresh = 1000,
  init = 0
)
print(fit_ref$time()$total)
