
functions {
#include functions/SIR.stan
#include functions/posterior.stan
}

data {
  int<lower=1> T;
  int y[T,1];
  real t0;
  real ts[T];
  real<lower=0> abs_tol_INF_;
  real<lower=0> rel_tol_INF_;
  int<lower=1> max_iter_INF_; 
  real<lower=0> abs_tol_REF_;
  real<lower=0> rel_tol_REF_;
  int<lower=1> max_iter_REF_;
}

transformed data {
  int N = 763;
  real x_r[0];
  int x_i[1] = {N};
  real y0[3] = {N-1.0, 1.0, 0.0};
}

parameters{
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> phi_inv;
}

transformed parameters {
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;
  real theta[2] = {beta, gamma};
  
  // Solve ODE
  real y_hat[T,3] = integrate_ode_rk45(SIR, y0, t0, ts, theta, x_r, x_i,
      abs_tol_INF_, rel_tol_INF_, max_iter_INF_);
  log_prior_na += log_prior_noadjustment(beta, gamma, phi_inv);
  log_lik_na += log_likelihood_noadjustment(y, y_hat, phi_inv, T);
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
  real y_hat_REF_[T,3] = integrate_ode_rk45(SIR, y0, t0, ts, theta, x_r, x_i,
      abs_tol_REF_, rel_tol_REF_, max_iter_REF_);
  log_prior_na_REF_ += log_prior_noadjustment(beta, gamma, phi_inv);
  log_lik_na_REF_ += log_likelihood_noadjustment(y, y_hat_REF_, phi_inv, T);
}

