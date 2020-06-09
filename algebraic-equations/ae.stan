
functions {
  vector system(vector y,        // unknowns
                vector theta,    // parameters
                real[] x_r,      // data (real)
                int[] x_i) {     // data (integer)
    vector[2] z;
    z[1] = y[1] - theta[1];
    z[2] = y[1] * y[2] - theta[2];
    return z;
  }
}

data {
  int<lower=1> n;
  real x[n, 2];
  real<lower=0> sigma;
  real<lower=0> rtol;
  int<lower=1> max_iter;
}

transformed data {
  vector[2] mu_guess = [1.0, 1.0]';
  real x_r[0];
  int x_i[0];
}

parameters {
  vector<lower=1>[2] theta;
}

transformed parameters {
  vector[2] mu;
  mu = algebra_solver(system, mu_guess, theta, x_r, x_i);
}

model{
  target += normal_lpdf(theta[1] | 1, 5);    // prior for alpha
  for(i in 1:n){
    for(j in 1:2){
        target += normal_lpdf(x[i,j] | mu[j], sigma);     // likelihood 
    }
  }
}

