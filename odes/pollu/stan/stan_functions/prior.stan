
// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real k_1) {
  real log_prior = 0.0;
  log_prior += normal_lpdf(sigma | 0.0, 1.0);
  log_prior += lognormal_lpdf(k_1 | 0.1, 1.0);
  return(log_prior);
}
