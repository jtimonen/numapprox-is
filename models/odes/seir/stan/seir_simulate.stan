
functions {
#include stan_functions/SEIR.stan
#include stan_functions/extract.stan
}

data {
  
  int<lower=1> K;         # number of age groups
  int<lower=1> D;         # number of days
  int<lower=1> S;         # number of ?
  int<lower=1> G;         # number of ?
  vector[K] age_dist;     # age distribution
  int<lower=1> pop_t;     # total population
  real<lower=0> ts[D];
  
  # Parameter values
  real<lower=0,upper=1> psi;       # proportion of symptomatic - Beta(58.3, 13.7)
  real<lower=0,upper=1> pii;       # proportion of exposed at t0 - Beta(1, 999)
  real<lower=0> phi[2];            # overdispersion for cases and deaths - Exp(0.01)
  vector[K]<lower=0,upper=1> eps;  # prop. of deaths among symptom. cases - Beta(1, 1)
  vector[K]<lower=0,upper=1> rho;  # reporting rate of symptom. cases - Beta(1, 1)
  
  real<lower=0> p_gamma[G];
  
  # BDF options
  real<lower=0> tol1;
  real<lower=0> tol2;
  int<lower=0> max_steps;
  
}

transformed data {
  real x_r[0];
  int x_i[0];
  real EPS = 1e-9;
  
  // Initial conditions
  for(k in 1:K){
    age_dist[k] = x_r[3+K*K + k];
    init[k] = age_dist[k] * (1-pii);
    init[K+k] = age_dist[k] * pii;
    init[2*K+k] = 0.0;
    init[3*K+k] = 0.0;
  }
}

model {
}

generated quantities {
  real y_hat[T,D] = integrate_ode_bdf(SEIR, y0, t0, ts, theta, x_r, x_i,
      tol1, tol2, max_steps);
      
  # Data to be sampled
  vector[D] out_icases;   // overall case incidence by day
  vector[D] out_ideaths;  // overall mortality incidence by day 
  vector[K] out_acases;   // final age distribution of cases (simplex)
  vector[K] out_adeaths;  // final age distribution of deaths (simplex)
  
  // Extract compartments from y
  vector[K] comp_S[S] = get_comp_S(y, age_dist, pii, pop_t, EPS);
  vector[K] comp_E[S] = get_comp_E(y, age_dist, pii, pop_t, EPS);
  vector[K] comp_I[S] = get_comp_I(y, pop_t, EPS, K);
  vector[K] comp_C[S+G] = get_comp_C(y, pop_t, EPS, K, G);
  vector[K] comp_diffC[S+G] = get_comp_diffC(comp_C, pop_t, EPS, S, G);
  vector[K] comp_diffM[S+G] = get_comp_diffM(comp_diffC, epsilon, p_gamma, EPS, S, G);
  vector[K] comp_M[S+G] = get_comp_M(comp_diffM, S, G);
  
  // Compute outcomes 
  for(i in t_data:S){
    out_icases[i-t_data+1] = sum(comp_diffC[i].*rho);
    out_ideaths[i-t_data+1] = sum(comp_diffM[i]);
  }
  out_acases = (comp_C[S,].*rho) ./ sum(comp_C[S,].*rho);
  out_adeaths = (comp_M[S,]) ./ sum(comp_M[S,]);

  // Evaluate likelihood
  for(i in 1:D) {
    icases[i]  = neg_binomial_2_rng(out_icases[i], out_icases[i] / phi[1] );
    ideaths[i] = neg_binomial_2_rng(out_ideaths[i], out_ideaths[i] / phi[2] );
  }
  acases  = multinomial_rng(out_acases);
  adeaths = multinomial_rng(out_adeaths);
  
}
