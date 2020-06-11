
# Setup
library(rstan)
library(loo)
source("functions.R")
model <- stan_model(file = 'sho_rk4.stan')
#model <- stan_model(file = 'sho_rk45.stan')
stan_data <- readRDS(file = 'simdat.rds')

# Additional data for reference method ode_integrate_rk45
stan_data$abs_tol_REF_  <- 1.0E-10
stan_data$rel_tol_REF_  <- 1.0E-10
stan_data$max_iter_REF_ <- 1.0E6

# Additional data for the fixed-step solver
step_size <- 0.5
stan_data <- add_interpolation_data(stan_data, step_size)

# Run sampling
fit <- sampling(object  = model,
                data    = stan_data,
                iter    = 2000,
                chains  = 4)

# Extract log posterior values (not Jacobian adjusted)
lh1 <- get_samples(fit, 'log_lik_na')
lh2 <- get_samples(fit, 'log_lik_na_REF_')
pr1 <- get_samples(fit, 'log_prior_na')
pr2 <- get_samples(fit, 'log_prior_na_REF_')
post1 <- lh1 + pr1
post2 <- lh2 + pr2

# PSIS
out <- psis(post1 - post2)
pareto_k <- out$diagnostics$pareto_k
cat('PARETO_K =', pareto_k, '\n')

# Get runtime
runtimes <-get_elapsed_time(fit)
print(runtimes)
