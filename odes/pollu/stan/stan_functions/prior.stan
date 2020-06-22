
// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real[] theta) {
  real log_prior = 0.0;
  log_prior += normal_lpdf(sigma | 0.0, 0.1);
  return(log_prior);
}
