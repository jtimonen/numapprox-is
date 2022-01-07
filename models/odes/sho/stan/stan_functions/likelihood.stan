// Likelihood without Jacobian adjustment
real log_likelihood_noadjustment(real[,] y, real[,] y_hat, real sigma, int T) {
  real log_lik = 0.0;
  for (t in 1:T){
    log_lik += normal_lpdf(y[t,1] | y_hat[t,1], sigma);
    log_lik += normal_lpdf(y[t,2] | y_hat[t,2], sigma);
  }
  return(log_lik);
}
