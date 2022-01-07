// Prior without Jacobian adjustment
real log_prior_noadjustment(real beta, real gamma, real phi_inv) {
  real log_prior = 0.0;
  log_prior += normal_lpdf(beta | 0.7, 0.5);
  log_prior += normal_lpdf(gamma | 0.2, 0.5);
  log_prior += exponential_lpdf(phi_inv | 0.2);
  return(log_prior);
}

// Likelihood without Jacobian adjustment
real log_likelihood_noadjustment(int[,] y, real[,] y_hat, real phi_inv, int T) {
  real log_lik = 0.0;
  for(i in 1:T){
      log_lik += neg_binomial_2_log_lpmf(y[i,1] | y_hat[i,2], phi_inv);
  }
  return(log_lik);
}
