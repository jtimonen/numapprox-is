  int<lower=1> T;
  real y[T,20];
  real t0;
  real ts[T];
  real<lower=0> k_data[25];
  
  real<lower=0> abs_tol_REF_;  // absolute tolerance for reference method
  real<lower=0> rel_tol_REF_;  // relative tolerance for reference method
  int<lower=1> max_iter_REF_;  // max number of iterations for reference method
  