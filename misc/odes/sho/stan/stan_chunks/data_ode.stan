  real<lower=0> abs_tol_REF_;  // absolute tolerance for BDF
  real<lower=0> rel_tol_REF_;  // relative tolerance for BDF
  int<lower=1> max_iter_REF_;  // max number of iterations for BDF
  
  real<lower=0> STEP_SIZE;
  int<lower=0> INTERP_R[T];
  real<lower=0> INTERP_A[T];
