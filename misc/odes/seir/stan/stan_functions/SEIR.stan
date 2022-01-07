// Sigmoidal switch
real switch_eta(real t, real t1, real eta, real nu, real xi) {
  return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
}

// ODE system
real[] SEIR(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
  int K = x_i[1];
  real tswitch = x_r[1];
  real dydt[(4*K)];  // SEIAR, then C and D
  real nI;           // total infectious
  real ntot;
  real beta;         // transmission rate
  real eta;          // reduction in transmission rate after quarantine
  real xi;           // slope of quarantine implementation
  real nu;           // shift of quarantine implementation
  real tau;          // incubation period
  real mu;           // infectious period
  real psi;          // probability of symptoms
  real p_tswitch;
  real contact[K*K]; // contact matrix, first K values
  // corresponds to number of contacts between 
  // age class 1 and other classes, etc
  
  real n_by_age[K];
  real f_inf[K];     // force of infection
  real init[K*4];
  real age_dist[K];
  real pii;           // number of cases at t0
  
  // Estimated parameters
  beta = theta[1];
  eta = theta[2];
  xi = theta[3];
  nu = theta[4];
  pii = theta[5];
  psi = theta[6];
  
  // Initial conditions
  for(k in 1:K){
    age_dist[k] = x_r[3+K*K + k];
    init[k] = age_dist[k] * (1-pii);
    init[K+k] = age_dist[k] * pii;
    init[2*K+k] = 0.0;
    init[3*K+k] = 0.0;
  }
  
  // Fixed parameters
  tau = 1.0/x_r[2];
  mu = 1.0/x_r[3];
  contact = x_r[4:(3+K*K)];
  
  // Total number of infectious people
  p_tswitch = switch_eta(t, tswitch, eta, nu, xi);
  
  // Force of infection by age classes: 
  // beta * p_tswitch * sum((#infected) / (#people) * (#contacts))
  for(k in 1:K){
    f_inf[k] = beta * p_tswitch * sum(to_vector(y[(2*K+1):(3*K)]) ./ to_vector(age_dist) .* to_vector(contact[(K*(k-1)+1):(k*K)]));
  }
  
  // Final dydt
  for (k in 1:K) {
    dydt[k]     = - f_inf[k] * (y[k]+init[k]);                                  // S
    dydt[K+k]   = f_inf[k] * (y[k]+init[k])- tau * (y[K+k]+init[K+k]);          // E
    dydt[2*K+k] = psi * tau * (y[K+k]+init[K+k]) - mu * (y[2*K+k]+init[2*K+k]); // I
    dydt[3*K+k] = psi * tau * (y[K+k]+init[K+k]);                               // C
  }
  return(dydt);
}

// Vector version of SEIR
vector odefun(real t, vector y, real[] theta, data real[] x_r, data int[] x_i){
  return to_vector(SEIR(t, to_array_1d(y), theta, x_r, x_i));
}
  