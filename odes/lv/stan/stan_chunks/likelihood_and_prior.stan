  // Compute prior and likelihood
  log_prior_na += log_prior_noadjustment(sigma, theta);
  log_lik_na += log_likelihood_noadjustment(y, y_hat, sigma, T);