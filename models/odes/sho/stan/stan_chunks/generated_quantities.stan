  real log_prior_na_REF_ = 0.0;
  real log_lik_na_REF_ = 0.0;
  
  // Solve ODE using reference method
  real y_hat_REF_[T,2] = integrate_ode_rk45(sho, y0, t0, ts, theta, x_r, x_i,
      abs_tol_REF_, rel_tol_REF_, max_iter_REF_);
  
  // Compute prior and likelihood
  log_prior_na_REF_ += log_prior_noadjustment(sigma, theta[1], y0);
  log_lik_na_REF_ += log_likelihood_noadjustment(y, y_hat_REF_, sigma, T);