// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real[] theta) {
  real log_prior = 0.0;
  log_prior += normal_lpdf(sigma | 1, 1);
  log_prior += gamma_lpdf(theta[1] | 4, 1);
  log_prior += gamma_lpdf(theta[2] | 4, 1);
  return(log_prior);
}

// Likelihood without Jacobian adjustment
real log_likelihood_noadjustment(real[,] y, real[,] y_hat, real sigma, int T) {
  real log_lik = 0.0;
  for (t in 1:T){
    log_lik += normal_lpdf(y[t,1] | y_hat[t,1], sigma);
    log_lik += normal_lpdf(y[t,2] | y_hat[t,2], sigma);
  }
  return(log_lik);
}
