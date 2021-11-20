// gsir_post.stan
functions {
  // Age-structured SIR system right-hand side
  vector SIR(real t, vector y, real beta, vector gamma, data matrix contacts,
      data vector pop_sizes) {
    int G = num_elements(pop_sizes);
    vector[2*G] dy_dt; // first G are susceptible, next G are infected
    vector[G] infection_rates;
    vector[G] recovery_rates;
    vector[G] lambda = rep_vector(0.0, G);
    for(g in 1:G){
      for(h in 1:G) {
        lambda[g] += contacts[g,h] * y[G+h]/pop_sizes[h];
      }
    }
    for(g in 1:G){
      dy_dt[g] = -  beta * lambda[g] * y[g];
      dy_dt[G+g] =  beta * lambda[g] * y[g] - gamma[g] * y[G+g];
    }
    return dy_dt;
  }
}

data {
  int<lower=1> N;                 // number of time points
  real t[N];                      // time points
  int<lower=1> G;                 // number of groups
  vector[G] pop_sizes;            // population sizes in each group
  vector[G] I0;                   // initial number of infected in each group
  matrix[G, G] contacts;          // contact matrix
  real<lower=0> rel_tol;          // ODE solver relative tolerance
  real<lower=0> abs_tol;          // ODE solver absolute tolerance
  int<lower=0> max_num_steps;     // ODE solver maximum number of steps
  int<lower=0> I_data[N, G];
}

transformed data {
  real t0 = 0.0;
  vector[2*G] x0;
  for(g in 1:G){
    x0[g] = pop_sizes[g] - I0[g]; // S
  }
  for(g in 1:G){
    x0[G + g] = I0[g]; // I
  }
}

parameters {
  real<lower=0> beta;
  vector<lower=0>[G] gamma;
  real<lower=0> phi_inv;
}

transformed parameters {
  real phi = inv(phi_inv);
}

model {
  beta ~ normal(2, 1);
  gamma ~ normal(0.3, 0.3);
  phi_inv ~ exponential(5);
  vector[2*G] x[N] = ode_rk45_tol(SIR, x0, t0, t, rel_tol, abs_tol, 
    max_num_steps, beta, gamma, contacts, pop_sizes);

  for(n in 1:N) {
    for(g in 1:G) {
      I_data[n,g] ~ neg_binomial_2(x[n][G+g] + 10*abs_tol, phi);
    }
  }
  
}
