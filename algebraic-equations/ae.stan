
functions {
  vector system(vector y, vector theta, real[] x_r, int[] x_i) {
    return [-y[1]^5 + y[1] + theta[1]]';
  }
}

data {
  int<lower=1> n;
  vector[n] x;
  real<lower=0> sigma;
  real<lower=0> rel_tol;
  int<lower=1> max_steps;
  vector[1] mu_guess;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  vector<lower=1>[1] theta;
}

transformed parameters {
  vector[1] mu;
  mu = algebra_solver(system, mu_guess, theta, x_r, x_i, rel_tol, 1e-6, max_steps);
}

model{
  target += normal_lpdf(theta[1] | 0, 5);    // prior
  target += normal_lpdf(x | mu[1], sigma);   // likelihood 
}

