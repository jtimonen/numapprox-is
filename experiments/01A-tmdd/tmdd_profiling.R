# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
source("../R/functions.R")
source("tmdd_setup.R")

cmdstanr::set_cmdstan_path("/Users/juhotimonen/Work/research/stan/cmdstan")
ITER <- 200
CHAINS <- 2
model <- tmdd_model(profile_oderhs = TRUE)

#  Helper function
get_profile_matrix <- function(fit) {
  p <- sapply(fit$cmdstanr_fit$profiles(), function(x) x)
  as.matrix(p[3:nrow(p), ])
}

# tol = 0.05
solver1 <- bdf(tol = 0.05, max_num_steps = 1e5)
fit1 <- model$sample(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solver = solver1, step_size = 0.1,
  iter_warmup = ITER,
  iter_sampling = ITER,
  chains = CHAINS
)
profs1 <- get_profile_matrix(fit1)

# tol = 0.02
solver2 <- bdf(tol = 0.02, max_num_steps = 1e5)
fit2 <- model$sample(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solver = solver2, step_size = 0.1,
  iter_warmup = ITER,
  iter_sampling = ITER,
  chains = CHAINS
)
profs2 <- get_profile_matrix(fit2)

# autodiff calls
adc1 <- c(24771063, 37775529, 34026450, 35440701)
adc2 <- c(14514573, 14042580, 16244124, 16331904)
