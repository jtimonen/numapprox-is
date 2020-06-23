// Air pollution model from: https://archimede.dm.uniba.it/~testset/problems/pollu.php

functions {
#include stan_functions/POLLU.stan
}

data {
  int<lower=1> T;
  real t0;
  real ts[T];
  real theta[25];
  real<lower=0> sigma;
  real<lower=0> tol1;
  real<lower=0> tol2;
  int<lower=0> max_steps;
}

transformed data {
  int D = 20;
  real x_r[0];
  int x_i[0];
  real y0[D] = {0, 0.2, 0, 0.04, 0, 0, 0.1, 0.3, 0.01, 0,
                0, 0, 0, 0, 0, 0, 0.007, 0, 0, 0};
}

model {
}

generated quantities {
  real y_hat[T,D] = integrate_ode_bdf(POLLU, y0, t0, ts, theta, x_r, x_i,
      tol1, tol2, max_steps);
  real y[T,D];
  // add measurement error
  for (t in 1:T) {
    for(d in 1:D){
      y[t, d] = y_hat[t, d] * (1.0 + normal_rng(0, sigma));
    }
  }
}
