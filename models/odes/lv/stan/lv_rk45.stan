
functions {
#include functions/LV.stan
#include functions/posterior.stan
}

data {
  int<lower=1> T;
  real y[T,2];
  real t0;
  real ts[T];
  real<lower=0> abs_tol_REF_;
  real<lower=0> rel_tol_REF_;
  int<lower=1> max_iter_REF_;
  
  real<lower=0> abs_tol_INF_;
  real<lower=0> rel_tol_INF_;
  int<lower=1> max_iter_INF_; 
}

transformed data {
  real x_r[0];
  int x_i[0];
  real y0[2] = {1.0, 1.0};
}

parameters{
  real<lower=0> sigma;
  real<lower=0> theta[2];
}

transformed parameters {
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;

  // Solve ODE
  real y_hat[T,2] = integrate_ode_rk45(LV, y0, t0, ts, theta, x_r, x_i,
      abs_tol_INF_, rel_tol_INF_, max_iter_INF_);
  log_prior_na += log_prior_noadjustment(sigma, theta);
  log_lik_na += log_likelihood_noadjustment(y, y_hat, sigma, T);
}

model {
  // Evaluate lp__ with the (invisible) Jacobian adjustment term included
  target += log_prior_na;
  target += log_lik_na;
}

generated quantities{
  real log_prior_na_REF_ = 0.0;
  real log_lik_na_REF_ = 0.0;
  // Solve ODE using reference method
  real y_hat_REF_[T,2] = integrate_ode_rk45(LV, y0, t0, ts, theta, x_r, x_i,
      abs_tol_REF_, rel_tol_REF_, max_iter_REF_);
  log_prior_na_REF_ += log_prior_noadjustment(sigma, theta);
  log_lik_na_REF_ += log_likelihood_noadjustment(y, y_hat_REF_, sigma, T);
}

