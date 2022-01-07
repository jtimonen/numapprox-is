// SIR infectious disease model

functions {
#include functions/SIR.stan
}

data {
  int<lower=1> T;
  real t0;
  real ts[T];
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
  real<lower=0> abs_tol;
  real<lower=0> rel_tol;
  int<lower=0> max_iter;
}

transformed data {
  int N = 763;
  real x_r[0];
  int x_i[1] = {N};
  real y0[3] = {N-1.0, 1.0, 0.0};
  real theta[2] = {beta, gamma};
}

model {
}

generated quantities {
  real y_hat[T,3] = integrate_ode_rk45(SIR, y0, t0, ts, theta, x_r, x_i,
      abs_tol, rel_tol, max_iter);
  int y[T] = neg_binomial_2_rng(col(to_matrix(y_hat), 2) + 1e-8, 1/phi_inv);
}
