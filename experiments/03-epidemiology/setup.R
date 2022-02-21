#!/usr/bin/env Rscript
# Setup SIRD/SEIR experiment

# R functions and requirements
source("../R/utils.R")
source("../R/functions.R")
library(odemodeling)


# Swiss data
load_data_switzerland <- function(parent_dir = ".") {
  fp <- file.path(parent_dir, "switzerland.rds")
  readRDS(file = fp)
}

# Create the model
seir_model <- function(prior_only = FALSE, ...) {
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
    real p_infected_survey = 0.0; // proportion of people having been infected at week 5
    for (i in 1:n_days-1){
      incidence[i] = -(y_sol[i+1][2] - y_sol[i][2] + y_sol[i+1][1] - y_sol[i][1]) * p_reported + delta;
    }
    for (i in t_survey_start:t_survey_end) {
      p_infected_survey += y_sol[i][4];
    }
    p_infected_survey = p_infected_survey / ((t_survey_end - t_survey_start + 1)*(pop_size));
    log_lik += binomial_lpmf(n_infected_survey | n_tested_survey, p_infected_survey + 1e-6);
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


# Compile model
model <- seir_model()

# Load data
dat <- load_data_switzerland("../../data/switzerland/")

# Create actual Stan data
add_data <- dat[c(
  "n_days", "cases", "tswitch", "t_survey_start",
  "t_survey_end", "n_infected_survey", "n_tested_survey"
)]

add_data$pop_size <- dat$N
add_data$D <- 4
add_data$n_days_m1 <- add_data$n_days - 1
t <- dat$ts
t0 <- dat$t0

# Other
init <- 0
step_size <- 0.1
