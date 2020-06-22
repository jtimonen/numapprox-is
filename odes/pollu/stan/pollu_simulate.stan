// Air pollution model from [1].
//  
// [1] Ernst Hairer and Gerhard Wanner.
//     Solving Ordinary Differential Equations II - 
//     Stiff and Differential-Algebraic Problems. Springer, 1991.
//

functions {
#include stan_functions/pollu.stan
}

data {
  int<lower=1> T;
  real t0;
  real ts[T];
  real theta[25];
  real<lower=0> sigma;
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
  real y_hat[T,D] = integrate_ode_bdf(POLLU, y0, t0, ts, theta, x_r, x_i);
  real y[T,D];
  // add measurement error
  for (t in 1:T) {
    for(d in 1:D){
      y[t, d] = y_hat[t, d] + normal_rng(0, sigma);
    }
  }
}
