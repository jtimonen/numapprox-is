  // These store non-jacobian-adjusted versions of prior and likelihood
  real log_prior_na = 0.0;
  real log_lik_na = 0.0;
  
  // This stores the ODE solution
  real y_hat[T,2];