// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real[] theta) {
  real log_prior = 0.0;
  log_prior += lognormal_lpdf(sigma | 1, 1);
  log_prior += lognormal_lpdf(theta[1] | 1, 1);
  log_prior += lognormal_lpdf(theta[2] | 1, 1);
  return(log_prior);
}

// Likelihood without Jacobian adjustment
real log_likelihood_noadjustment(real[,] y, real[,] y_hat, real sigma, int T) {
  real log_lik = 0.0;
  log_lik += normal_lpdf(y[:, 1] | y_hat[:, 1], sigma);
  log_lik += normal_lpdf(y[:, 2] | y_hat[:, 2], sigma);
  return(log_lik);
}
