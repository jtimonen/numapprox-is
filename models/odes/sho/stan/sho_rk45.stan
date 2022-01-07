
functions {
#include stan_functions/sho.stan
#include stan_functions/likelihood.stan
#include stan_functions/prior.stan
}

data {
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
  y_hat = integrate_ode_rk45(sho, y0, t0, ts, theta, x_r, x_i);
#include stan_chunks/likelihood_and_prior.stan
}

model {
#include stan_chunks/model.stan
}

generated quantities{
#include stan_chunks/generated_quantities.stan
}

