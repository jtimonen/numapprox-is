
// Prior without Jacobian adjustment
real log_prior_noadjustment(real sigma, real theta, real[] y0) {
  real log_prior = 0.0;
  log_prior += cauchy_lpdf(sigma | 0.0, 2.5);
  log_prior += normal_lpdf(theta | 0.0, 1.0);
  log_prior += normal_lpdf(y0[1]  | 0.0, 1.0);
  log_prior += normal_lpdf(y0[2]  | 0.0, 1.0);
  return(log_prior);
}
