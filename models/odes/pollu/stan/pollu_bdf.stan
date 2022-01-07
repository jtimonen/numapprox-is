
functions {
#include stan_functions/POLLU.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan
}

data {
  real<lower=0> abs_tol_INF_;  // absolute tolerance for reference method
  real<lower=0> rel_tol_INF_;  // relative tolerance for reference method
  int<lower=1> max_iter_INF_;  // max number of iterations for reference method
#include stan_chunks/data.stan
#include stan_chunks/data_ode.stan
}

transformed data {
#include stan_chunks/transformed_data.stan
}

parameters{
#include stan_chunks/parameters.stan
}

transformed parameters {
#include stan_chunks/transformed_parameters.stan
  y_hat = integrate_ode_bdf(POLLU, y0, t0, ts, theta, x_r, x_i,
      abs_tol_INF_, rel_tol_INF_, max_iter_INF_);
#include stan_chunks/likelihood_and_prior.stan
}

model {
#include stan_chunks/model.stan
}

generated quantities{
#include stan_chunks/generated_quantities.stan
}

