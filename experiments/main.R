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

# Options
stan_opts <- list(
  sig_figs = 18 # number of significant figures to store in floats
)

# R functions
source("functions.R")
source("setup_sir.R")

# Create Stan models
dat <- setup_standata_sir()
code_parts <- setup_stancode_sir(solver = "rk45")
stanmodels <- create_cmdstan_models(code_parts)

# Fit prior model
prior_fit <- stanmodels$prior$sample(
  sig_figs = stan_opts$sig_figs,
  seed = stan_opts$seed,
  refresh = 1000
)
prior_draws <- prior_fit$draws(c("beta", "gamma", "phi_inv"))

# Simulate solutions using prior draws
solver_args_gen <- list(
  rel_tol = 1e-12,
  abs_tol = 1e-12,
  max_num_steps = 1e9
)
prior_sim <- simulate(
  stanmodels$simulator, prior_draws, dat,
  solver_args_gen, stan_opts
)
y_gen_rvar <- posterior::as_draws_rvars(prior_sim$draws("y_gen"))$y_gen
y_gen_arr <- posterior::as_draws_array(y_gen_rvar)

df <- data.frame(dat$t, t(y_gen_rvar))
colnames(df) <- c("Day", "S", "I")

plt_S <- ggplot(df, aes(x = Day, y = S)) +
  stat_lineribbon(alpha = 0.6) +
  scale_fill_brewer() +
  ggtitle("Number of susceptible (S), population size = 763")
plt_I <- ggplot(df, aes(x = Day, y = I)) +
  stat_lineribbon(alpha = 0.6) +
  scale_fill_brewer() +
  ggtitle("Number of infected (I), population size = 763")

# Plot generated data
i_draws <- as_draws_array(y_gen_rvar[2, ])
NS <- dim(i_draws)[1]
plot(dat$t, i_draws[1, 1, ], ylim = c(0, 2000), ylab = "Infected", xlab = "Day",
     pch=".")
for (s in 1:NS) {
  points(dat$t, i_draws[s, 1, ], pch = "x", col = scales::alpha("black", 0.1))
}
lines(c(0, 14), c(dat$pop_size, dat$pop_size), col = "firebrick", lty = 2)


# Plot solutions
title <- "Solutions using prior param draws"
plot_sir_example_solutions(prior_sim, dat, thin = 10, main = title)

## Analyzing the numerical method
TOLS <- 10^seq(-9, -4)
atols <- TOLS
rtols <- TOLS
sims <- simulate_many(
  stanmodels$sim, prior_draws, dat, stan_opts,
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
fit <- sample_posterior(stanmodels$posterior, dat,
  solver_args_sample, stan_opts,
  refresh = 1000,
  init = 0
)

print(fit$time())
print(fit)
post_draws <- fit$draws(c("beta", "gamma", "phi_inv"))

# Simulate and plot generated solutions
post_sim <- simulate(stanmodels$sim, post_draws, dat, solver_args_sample, stan_opts)

title <- "Solutions using posterior param draws"
plot_sir_example_solutions(post_sim, dat, thin = 10, main = title)

# Tune the numerical method in $M_{high}$ so that it is reliable at these draws.
TOLS <- 10^seq(-4, -12) # could be just halving for example
tuning_tol <- 0.0001
tuning <- tune_solver(
  TOLS, stanmodels$sim, post_draws, dat, stan_opts,
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
fit_ref <- sample_posterior(stanmodels$posterior, dat,
  solver_args_refsample, stan_opts,
  refresh = 1000,
  init = 0
)
print(fit_ref$time()$total)
