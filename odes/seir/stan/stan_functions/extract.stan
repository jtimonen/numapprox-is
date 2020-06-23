// Extract compartment S
vector[] get_comp_S(real[,] y, vector age_dist, real pii, data int pop_t, data real EPS){
  int S = size(y);
  int K = num_elements(age_dist);
  vector[K] comp[S];
  for(i in 1:S) {
    comp[i] = (to_vector(y[i,1:K]) + to_vector(age_dist) * (1-pii) + EPS) * pop_t;
  }
  return(comp);
}

// Extract compartment E
vector[] get_comp_E(real[,] y, vector age_dist, real pii, data int pop_t, data real EPS){
  int S = size(y);
  int K = num_elements(age_dist);
  vector[K] comp[S];
  for(i in 1:S) {
    comp[i] = (to_vector(y[i,(K+1):(2*K)]) + to_vector(age_dist) * pii + EPS) * pop_t;
  }
  return(comp);
}

// Extract compartment I
vector[] get_comp_I(real[,] y, data int pop_t, data real EPS, int K){
  int S = size(y);
  vector[K] comp[S];
  for(i in 1:S) {
    comp[i] = (to_vector(y[i,(2*K+1):(3*K)]) + EPS) * pop_t;
  }
  return(comp);
}

// Extract compartment C
vector[] get_comp_C(real[,] y, data int pop_t, data real EPS, int K, int G){
  int S = size(y);
  vector[K] comp[S+G];
  for(i in 1:S) {
    comp[i] = (to_vector(y[i,(3*K+1):(4*K)]) + EPS) * pop_t;
  }
  for(g in 1:G){
    comp[S+g] = comp[S];
  }
  return(comp);
}

// Get diffC (lagged difference of cumulative incidence of symptomatics)
vector[] get_comp_diffC(vector[] comp_C, data int pop_t, data real EPS, int S, int G){
  int K = num_elements(comp_C[1]);
  vector[K] comp[S+G];
  for(i in 1:S) {
    comp[i] = i==1 ? comp_C[i,] : EPS*pop_t + comp_C[i,] - comp_C[i-1,]; 
  }
  for(g in 1:G){
    comp[S+g] = rep_vector(EPS, K);
  }
  return(comp);
}

// Get diffM (difference in mortality)
vector[] get_comp_diffM(vector[] comp_diffC, vector epsilon, data real[] p_gamma, data real EPS, int S, int G){
  int K = num_elements(comp_diffC[1]);
  vector[K] comp[S+G];
  for(i in 1:(G+S)){
    comp[i] = rep_vector(EPS, K);
  }
  for(i in 1:S) {
    for(g in 1:G) {
      comp[i+g] += comp_diffC[i] .* epsilon * p_gamma[g] ;
    }
  }
  return(comp);
}

// Get comp_M (cumulative mortality)
vector[] get_comp_M(vector[] comp_diffM, int S, int G){
  int K = num_elements(comp_diffM[1]);
  vector[K] comp[S+G];
  for(i in 1:(S+G)) {
    for(k in 1:K) {
      comp[i,k] = sum(comp_diffM[1:i,k]);
    }
  }
  return(comp);
}

// Get comp_D
vector[] get_comp_D(vector[] comp_C, real psi, int K){
  int S = size(comp_C);
  vector[K] comp[S];
  for(i in 1:S){
    comp[i] = (1.0-psi)/psi * comp_C[i];
  }
  return(comp);
}

// Get comp_diffD
vector[] get_comp_diffD(vector[] comp_diffC, real psi, int K){
  int S = size(comp_diffC);
  vector[K] comp[S];
  for(i in 1:S){
    comp[i] = (1.0-psi)/psi * comp_diffC[i];
  }
  return(comp);
}


