# Data
setup_data_sir <- function(){
  I_data <- outbreaks::influenza_england_1978_school$in_bed
  N <- length(I_data)
  data_list <- list(
    I_data = I_data,
    N = N,
    t = seq(1, N),
    I0 = 1,
    pop_size = 763,
    D = 3
  )
  return(data_list)
}

setup_stancode_sir <- function(solver = "rk45"){
  
  pars <- "
    real<lower=0> beta;
    real<lower=0> gamma;
    real<lower=0> phi_inv;
  "
  tpars <- "  real phi = inv(phi_inv);"
  prior <- "
    beta ~ normal(2, 1);
    gamma ~ normal(0.4, 0.5);
    phi_inv ~ exponential(5);
  "
  funs <- "
  // SIR system right-hand side
  vector SIR(real t, vector y, data int[] a0, vector theta) {
    vector[3] dx_dt;
    int pop_size = a0[1];
    real infection_rate = theta[1] * y[2] * y[1] / pop_size;
    real recovery_rate = theta[2] * y[2];
    dx_dt[1] = - infection_rate;
    dx_dt[2] = infection_rate - recovery_rate;
    dx_dt[3] = recovery_rate;
    return dx_dt;
  }
  "
  data <- "
    int<lower=1> N; // number of time points
    int<lower=3,upper=3> D; // dimension of ODE
    real t[N]; // time points
    int<lower=1> pop_size; // population size
    real<lower=1> I0; // initial number of infected
    real<lower=0> rel_tol; // ODE solver relative tolerance
    real<lower=0> abs_tol; // ODE solver absolute tolerance
    int<lower=0> max_num_steps; // ODE solver maximum number of steps
  "
  tdata <- "
    real t0 = 0.0;
    int a0[1] = {pop_size};
    vector[D] x0 = to_vector({pop_size - I0, I0, 0.0}); // {S, I, R}
  "
  obsdata <- "  int<lower=0> I_data[N];"
  if(solver=="rk45") {
    odesolve <- "  vector[D] x[N] = ode_rk45_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, a0, to_vector({beta, gamma}));"
  } else if(solver=="bdf") {
    odesolve <- "  vector[D] x[N] = ode_bdf_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, a0, to_vector({beta, gamma}));"
  } else {
    stop("Invalid solver argument!")
  }

  likelihood <- "
    for(n in 1:N) {
      I_data[n] ~ neg_binomial_2(x[n][2] + 10*abs_tol, phi);
    }
  "
  genquant <- "
    int y[D, N];
    real log_lik = 0.0;
    for(n in 1:N) {
      for(d in 1:D) {
        y[d,n] = neg_binomial_2_rng(x[n][d] + 10*abs_tol, phi); 
      }
      log_lik += neg_binomial_2_lpmf(I_data[n] | x[n][2] + 10*abs_tol, phi);
    }
  "
  
  # Return
  list(
    functions = funs,
    pars = pars,
    tpars = tpars,
    prior = prior,
    data = data,
    tdata = tdata,
    obsdata = obsdata,
    odesolve = odesolve,
    likelihood = likelihood,
    genquant = genquant
  )
}
