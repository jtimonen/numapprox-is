library(odemodeling)

# Lotka-Volterra ------------------------------------------------------------

# Create Lotka-Volterra model
ode_model_lv <- function(prior_only = FALSE, ...) {

  # Dimensions and other data
  N <- stan_dim("N", lower = 1) # number of time points
  D <- stan_dim("D", lower = 1) # ODE system dimension
  y_obs <- stan_vector_array("y_obs", length = D, dims = list(N))
  y_obs_init <- stan_vector("y_obs_init", length = D)

  # Define ODE system parameters and their priors
  lv_par <- list(
    stan_param(stan_var("alpha", lower = 0), "normal(1, 0.5)"),
    stan_param(stan_var("beta", lower = 0), "normal(0.05, 0.05)"),
    stan_param(stan_var("gamma", lower = 0), "normal(1, 0.5)"),
    stan_param(stan_var("delta", lower = 0), "normal(0.05, 0.05)")
  )

  # Define noise parameter and its prior
  sigma_par <- stan_param(
    stan_vector("sigma", lower = 0, length = D),
    "lognormal(-1, 1)"
  )

  # Define initial point as parameter and its prior
  y0_par <- stan_param(
    stan_vector("y0", length = D, lower = 0),
    "lognormal(log(10), 1)"
  )

  # Initial point on log scale
  log_y0 <- stan_transform(
    stan_vector("log_y0", length = D),
    "parameters",
    "log(y0)"
  )

  # Define ODE system right-hand side
  odefun_body <- "
    real u = y[1]; // predator
    real v = y[2]; // prey
    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return to_vector({du_dt, dv_dt});
  "

  # Define log-likelihood function body
  loglik_body <- "
    real loglik = lognormal_lpdf(y_obs_init | log_y0, sigma);
    for(n in 1:N) {
      loglik += lognormal_lpdf(y_obs[n] | log(y_sol[n]), sigma);
    }
    return(loglik);
  "

  # Set loglik depending on whether creating only prior model
  if (prior_only) {
    loglik_body <- ""
    loglik_vars <- list(log_y0, sigma_par)
  } else {
    loglik_vars <- list(log_y0, sigma_par, y_obs, y_obs_init)
  }

  # Return
  odemodeling::ode_model(
    N = N,
    odefun_vars = lv_par,
    odefun_body = odefun_body,
    odefun_init = y0_par,
    loglik_vars = loglik_vars,
    loglik_body = loglik_body,
    ...
  )
}

# TMDD --------------------------------------------------------------------

# Create target mediated drug disposition model
ode_model_tmdd <- function(prior_only = FALSE, ...) {

  # Dimensions and other data
  N <- stan_dim("N", lower = 1) # number of time points
  D <- stan_dim("D", lower = 1) # ODE system dimension
  L0 <- stan_var("L0", lower = 0) # initial bolus
  P_obs <- stan_vector("P_obs", length = N) # observations of P

  # Define kinetic parameters and their priors
  k_par <- list(
    stan_param(stan_var("k_on", lower = 0), "lognormal(-1, 0.3)"),
    stan_param(stan_var("k_off", lower = 0), "lognormal(0, 0.3)"),
    stan_param(stan_var("k_in", lower = 0), "lognormal(0, 0.3)"),
    stan_param(stan_var("k_out", lower = 0), "lognormal(0, 0.3)"),
    stan_param(stan_var("k_eL", lower = 0), "lognormal(-1, 0.3)"),
    stan_param(stan_var("k_eP", lower = 0), "lognormal(-3, 0.3)")
  )

  # Define noise parameter and its prior
  sigma_par <- stan_param(stan_var("sigma", lower = 0), "lognormal(1, 0.3)")

  # Define transformed parameters
  R0 <- stan_transform(stan_var("R0"), "parameters", "k_in/k_out")
  y0 <- stan_transform(
    decl = stan_vector("y0", length = D),
    origin = "parameters",
    code = "to_vector({L0, R0, 0.0})"
  )

  # Define ODE system right-hand side
  odefun_body <- "
    vector[3] dy_dt; // L, R, P
    real L = y[1];
    real R = y[2];
    real P = y[3];
    real rem = k_on*L*R - k_off*P;
    dy_dt[1] = - k_eL*L - rem;
    dy_dt[2] = k_in - k_out*R - rem;
    dy_dt[3] = rem - k_eP*P;
    return dy_dt;
  "

  # Define log-likelihood function body
  loglik_body <- "
    real loglik = 0.0;
    for(n in 1:N) {
      loglik += normal_lpdf(P_obs[n] | y_sol[n][3], sigma);
    }
    return(loglik);
  "

  # Set loglik depending on whether creating only prior model
  if (prior_only) {
    loglik_body <- ""
    loglik_vars <- list(sigma_par) # no P_obs
  } else {
    loglik_vars <- list(sigma_par, P_obs)
  }

  # Return
  odemodeling::ode_model(
    N = N,
    odefun_vars = k_par,
    odefun_body = odefun_body,
    odefun_init = y0,
    loglik_vars = loglik_vars,
    loglik_body = loglik_body,
    other_vars = list(L0, R0),
    ...
  )
}


# SEIR ---------------------------------------------------------------------

ode_model_seir <- function(prior_only = FALSE, ...) {
  if (prior_only) {
    stop("prior_only = FALSE not implemented!")
  }

  # Data block
  n_days <- stan_dim("n_days", lower = 1)
  n_days_m1 <- stan_dim("n_days_m1", lower = 0)
  tswitch <- stan_var("tswitch")
  pop_size <- stan_var("pop_size") # population size
  cases <- stan_array("cases", dims = list(n_days), type = "int")
  delta <- stan_transform(stan_var("delta"), "data", "10*abs_tol")

  # Antibody survey data
  t_survey_start <- stan_var("t_survey_start", type = "int")
  t_survey_end <- stan_var("t_survey_end", type = "int")
  n_infected_survey <- stan_var("n_infected_survey", type = "int")
  n_tested_survey <- stan_var("n_tested_survey", type = "int")

  # SEIR parameters
  gamma <- stan_param(stan_var("gamma", lower = 0), prior = "normal(0.4, 0.5)")
  beta <- stan_param(stan_var("beta", lower = 0), prior = "normal(2, 1)")
  a <- stan_param(stan_var("a", lower = 0), prior = "normal(0.4, 0.5)")

  # Observation model parameters phi_inv
  phi_inv_var <- stan_var("phi_inv", lower = 0)
  phi_var <- stan_var("phi", lower = 0)
  phi_inv <- stan_param(phi_inv_var, "exponential(5);")
  phi <- stan_transform(phi_var, "parameters", "inv(phi_inv);")

  # Initial infected and exposed
  i0 <- stan_param(stan_var("i0", lower = 0), prior = "normal(0, 10)")
  e0 <- stan_param(stan_var("e0", lower = 0), prior = "normal(0, 10)")

  # Reporting rate (probability for an infected person to be reported)
  p_reported <- stan_param(
    decl = stan_var("p_reported", lower = 0, upper = 1),
    prior = "beta(1, 2)"
  )

  # Slope of quarantine implementation
  xi_raw <- stan_param(
    decl = stan_var("xi_raw", lower = 0, upper = 1), prior = "beta(1, 1)"
  )

  # Reduction in transmission due to control measures
  eta <- stan_param(
    decl = stan_var("eta", lower = 0, upper = 1), prior = "beta(2.5, 4)"
  )

  # Shift of quarantine implementation
  nu <- stan_param(
    decl = stan_var("nu", lower = 0), prior = "exponential(1./5)"
  )

  # Transformed xi
  xi <- stan_transform(
    decl = stan_var("xi"), origin = "parameters", code = "xi_raw + 0.5"
  )

  # All ODE parameters
  ode_params <- list(gamma, beta, a, eta, nu, xi, tswitch)

  # Initial state
  D <- stan_dim("D", lower = 4, upper = 4) # SEIR
  y0 <- stan_transform(
    decl = stan_vector("y0", length = D),
    origin = "data",
    code = "to_vector({0.0, 0.0, 0.0, 0.0})"
  )

  # All odefun variables
  odefun_vars <- c(ode_params, list(pop_size, e0, i0))

  # All loglik variables
  loglik_vars <- list(
    cases, phi, p_reported, t_survey_start, t_survey_end,
    n_tested_survey, n_infected_survey, pop_size, delta
  )

  # Incidence
  incidence_gq <- stan_transform(
    stan_vector("incidence_gq", length = n_days_m1),
    origin = "model",
    code = "
      for (i in 1:n_days-1){
        incidence_gq[i] = -(y_sol_gq[i+1][2] - y_sol_gq[i][2] +
          y_sol_gq[i+1][1] - y_sol_gq[i][1]) * p_reported + delta;
      }"
  )

  # Other variables
  other_vars <- list(phi_inv, xi_raw, incidence_gq)

  # Function bodies
  odefun_body <- "
    real forcing = eta + (1 - eta) / (1 + exp(xi * (t - tswitch - nu)));
    real S = pop_size - e0 - i0 + y[1];
    real E = e0 + y[2];
    real I = i0 + y[3];
    real R = y[4];
    real exposing_rate = forcing * beta * I * S / pop_size;
    
    real dS_dt = - exposing_rate;
    real dE_dt = exposing_rate - a * E;
    real dI_dt = a * E - gamma * I;
    real dR_dt = gamma * I;
      
    return to_vector({dS_dt, dE_dt, dI_dt, dR_dt});
  "

  loglik_body <- "
    real log_lik = 0.0;
    real incidence[n_days - 1];
    real p_infected_survey; // proportion of people having been infected at week 5
    for (i in 1:n_days-1){
      incidence[i] = -(y_sol[i+1][2] - y_sol[i][2] + y_sol[i+1][1] - y_sol[i][1]) * p_reported + delta;
    }
    // mean number of infected + recovered people during week 5
    p_infected_survey = mean(to_vector(y_sol[t_survey_start:t_survey_end][4])) / pop_size;
    log_lik += binomial_lpmf(n_infected_survey | n_tested_survey, p_infected_survey);
    log_lik += neg_binomial_2_lpmf(cases[1:(n_days-1)] | incidence, phi);
    return(log_lik);
  "

  # Return
  ode_model(
    N = n_days,
    odefun_vars = odefun_vars,
    odefun_body = odefun_body,
    odefun_init = y0,
    loglik_vars = loglik_vars,
    loglik_body = loglik_body,
    other_vars = other_vars,
    ...
  )
}
