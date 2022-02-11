functions {
  // TMDD system right-hand side
  vector TMDD(real t, vector y, vector k) {
    vector[3] dy_dt; // L, R, P
    real L = y[1];
    real R = y[2];
    real P = y[3];
    real rem = k[1]*L*R - k[2]*P;
    dy_dt[1] = - k[5]*L - rem;
    dy_dt[2] = k[3] - k[4]*R - rem;
    dy_dt[3] = rem - k[6]*P;
    return dy_dt;
  }
  
}

data {
  int<lower=1> N;            // number of time points
  real t[N];                 // time points
  real L0;                   // initial bolus
  

  real<lower=0> rel_tol;          // ODE solver relative tolerance
  real<lower=0> abs_tol;          // ODE solver absolute tolerance
  int<lower=0> max_num_steps;     // ODE solver maximum number of steps

  vector[N] y;
}

transformed data {
    real t0 = 0.0;
  
}

parameters {
  vector<lower=0>[6] k; // on, off, in, out, eL, eP
  real<lower=0> sigma;
  
}

transformed parameters {
  
    real R0 = k[3]/k[4];
    vector[3] x0 = to_vector({L0, R0, 0.0});
  
}

generated quantities {
  vector[3] x[N] = ode_bdf_tol(TMDD, x0, t0, t, rel_tol, abs_tol, max_num_steps, k);

  vector[N] y_gen;
  real log_lik = 0.0;
  for(n in 1:N) {
    y_gen[n] = normal_rng(x[n][1], sigma);
    log_lik += normal_lpdf(y[n] | x[n][1], sigma);
  }
  
}
