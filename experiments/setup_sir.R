# Data
setup_standata_sir <- function(N = 14) {
  data_list <- list(
    N = N,
    t = as.numeric(seq(1, N)),
    I0 = 1,
    pop_size = 763,
    D = 2
  )
  return(data_list)
}

# Stan code parts
setup_stancode_sir <- function(solver = "rk45") {
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
    vector[2] dx_dt;
    int pop_size = a0[1];
    real infection_rate = theta[1] * y[2] * y[1] / pop_size;
    real recovery_rate = theta[2] * y[2];
    dx_dt[1] = - infection_rate;
    dx_dt[2] = infection_rate - recovery_rate;
    return dx_dt;
  }
  "
  data <- "
    int<lower=1> N; // number of time points
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
    vector[2] x0 = to_vector({pop_size - I0, I0}); // {S, I}
  "
  obsdata <- "  int<lower=0> I_data[N];"
  if (solver == "rk45") {
    odesolve <- "  vector[2] x[N] = ode_rk45_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, a0, to_vector({beta, gamma}));"
  } else if (solver == "bdf") {
    odesolve <- "  vector[2] x[N] = ode_bdf_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, a0, to_vector({beta, gamma}));"
  } else {
    stop("Invalid solver argument!")
  }

  likelihood <- "
    for(n in 1:N) {
      I_data[n] ~ neg_binomial_2(x[n][2] + 10*abs_tol, phi);
    }
  "
  genquant <- "
    int y_gen[2, N];
    for(n in 1:N) {
      y_gen[1][n] = neg_binomial_2_rng(x[n][1] + 10*abs_tol, phi);
      y_gen[2][n] = neg_binomial_2_rng(x[n][2] + 10*abs_tol, phi);
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

# Function that plots SIR solutions against data (of infected)
plot_sir_example_solutions <- function(sim, data, thin = 1, main = "") {
  N <- data$N
  x_sim <- posterior::thin_draws(posterior::merge_chains(sim$draws("x")), thin)
  I_sim <- x_sim[, 1, (N + 1):(2 * N), drop = TRUE]
  plot(data$t, rep(data$pop_size, data$N),
    type = "l", lty = 2,
    col = "gray70", ylim = c(0, 800), ylab = "Infected", xlab = "Day",
    main = main
  )
  for (i_draw in 1:nrow(I_sim)) {
    I <- as.vector(I_sim[i_draw, ])
    lines(data$t, I, col = scales::alpha("firebrick", 0.1))
  }
}
