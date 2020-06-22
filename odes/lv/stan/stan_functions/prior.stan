
// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real[] theta) {
  real log_prior = 0.0;
  log_prior += lognormal_lpdf(sigma | 0.1, 1.0);
  log_prior += lognormal_lpdf(theta[1] | 1.0, 1.0);
  log_prior += lognormal_lpdf(theta[2] | 1.0, 1.0);
  log_prior += lognormal_lpdf(theta[3] | 1.0, 1.0);
  return(log_prior);
}
