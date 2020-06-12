
# Setup
library(rstan)
library(loo)
source("functions.R")
#model <- stan_model(file = 'sho_rk4.stan')
model <- stan_model(file = 'sho_rk45.stan')

n_chains <- 4
SIGMA    <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
S        <- length(SIGMA)
PARETO_K <- rep(0, S)
RUNTIMES <- matrix(0, S, n_chains)

for(i in 1:S){
  
  # Load data
  sigma <- SIGMA[i]
  fn <- paste0('data/simdata_sigma_', sigma, '.rds')
  stan_data <- readRDS(file = fn)
  
  # Additional data for reference method ode_integrate_rk45
  stan_data$abs_tol_REF_  <- 1.0E-6
  stan_data$rel_tol_REF_  <- 1.0E-6
  stan_data$max_iter_REF_ <- 1.0E6
  
  # Additional data for the fixed-step solver
  step_size <- 0.5
  stan_data <- add_interpolation_data(stan_data, step_size)
  
  # Run sampling
  fit <- sampling(object  = model,
                  data    = stan_data,
                  iter    = 4000,
                  chains  = n_chains,
                  control = list(adapt_delta = 0.8))
  
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
  cat('sigma =', sigma, 'pareto_k =', pareto_k, '\n')
  PARETO_K[i] <- pareto_k
  
  # Get runtime
  runtimes <- get_elapsed_time(fit)
  RUNTIMES[i,] <- rowSums(runtimes)
  
}

print(PARETO_K)
