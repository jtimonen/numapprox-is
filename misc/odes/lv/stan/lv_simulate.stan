// Lotka-Volterra system

functions {
#include functions/LV.stan
}

data {
  int<lower=1> T;
  real y0[2];
  real t0;
  real ts[T];
  real theta[2];
  real<lower=0> sigma;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

model {
}

generated quantities {
  real y_hat[T,2] = integrate_ode_rk45(LV, y0, t0, ts, theta, x_r, x_i);
  real y[T,2];
  // add measurement error
  for (t in 1:T) {
    y[t, 1] = y_hat[t, 1] + normal_rng(0, sigma);
    y[t, 2] = y_hat[t, 2] + normal_rng(0, sigma);
  }
}
