// A phamacokinetic model

functions {
#include stan_functions/PKPD.stan
}

data {
  int<lower=1> T;
  real t0;
  real ts[T];
  real theta[5];
  real<lower=0> sigma;
}

transformed data {
  int D = 5;
  real x_r[0];
  int x_i[0];
  real y0[D] = {100.0, 0.0, 0.0, 0.0, 5.0};
}

model {
}

generated quantities {
  real y_hat[T,D] = integrate_ode_rk45(PKPD, y0, t0, ts, theta, x_r, x_i);
  real y[T,D];
  // add measurement error
  for (t in 1:T) {
    for(d in 1:D){
      y[t, d] = y_hat[t, d] + normal_rng(0, sigma);
    }
  }
}
