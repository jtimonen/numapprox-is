// Likelihood without Jacobian adjustment
real log_likelihood_noadjustment(real[,] y, real[,] y_hat, real sigma, 
    int T, int D) {
  real log_lik = 0.0;
  for (t in 1:T){
    for (d in 1:D){
      log_lik += normal_lpdf(y[t,d] | y_hat[t,d], sigma);
    }
  }
  return(log_lik);
}
