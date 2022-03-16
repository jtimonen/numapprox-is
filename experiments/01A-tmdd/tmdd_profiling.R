# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
source("../R/functions.R")
source("tmdd_setup.R")

ITER <- 2000
CHAINS <- 4
model <- tmdd_model(profile_oderhs = TRUE)

# tol = 0.05
solver1 <- bdf(tol = 0.05, max_num_steps = 1e5)
fit1 <- model$sample(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solver = solver1, step_size = 0.1,
  iter_warmup = ITER,
  iter_sampling = ITER,
  chains = CHAINS
)
profs1 <- sapply(fit1$cmdstanr_fit$profiles(), function(x) x)

# tol = 0.02
solver2 <- bdf(tol = 0.02, max_num_steps = 1e5)
fit2 <- model$sample(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solver = solver2, step_size = 0.1,
  iter_warmup = ITER,
  iter_sampling = ITER,
  chains = CHAINS
)
profs2 <- sapply(fit2$cmdstanr_fit$profiles(), function(x) x)

# autodiff calls
adc1 <- c(24771063, 37775529, 34026450, 35440701)
adc2 <- c(14514573, 14042580, 16244124, 16331904)
