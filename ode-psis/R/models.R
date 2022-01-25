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


# GSIR --------------------------------------------------------------------

ode_model_gsir <- function(prior_only = FALSE, ...) {

  # Time points
  N <- stan_dim("N", lower = 0) # number of timepoints

  # Data needed by ODE function
  G <- stan_dim("G", lower = 1) # number of groups
  pop_sizes <- stan_vector("pop_sizes", G) # population sizes in each group
  I0 <- stan_vector("I0", G, lower = 0) # initial no. infected in each group
  contacts <- stan_matrix("contacts", G, G) # contact matrix

  # ODE function parameters
  beta <- stan_param(stan_var("beta", lower = 0), "normal(2, 1)")
  gamma_decl <- stan_vector("gamma", lower = 0, length = G)
  mu_decl <- stan_vector("mu", lower = 0, length = G)
  gamma <- stan_param(gamma_decl, "normal(0.4, 0.5)")
  mu <- stan_param(mu_decl, "normal(0.4, 0.5)")

  # All odefun variables
  odefun_vars <- list(pop_sizes, I0, contacts, beta, gamma, mu)

  # Initial value
  D <- stan_dim("D", lower = 0) # number of ODE dimensions
  y0_var <- stan_vector("y0", length = D)
  y0_code <-
    "
    for(g in 1:G) {
      y0[g] = pop_sizes[g] - I0[g]; // initial number of S
    }
    for(g in 1:G) {
      y0[G + g] = I0[g]; // initial number of I
    }
    for(g in 1:G) {
      y0[2*G + g] = 0; // initial number of R
    }
    for(g in 1:G) {
      y0[3*G + g] = 0; // initial number of D
    }
  "
  y0 <- stan_transform(y0_var, "data", y0_code)

  # Observation model data
  delta <- stan_var("delta", lower = 0)
  deaths_cumulative <- stan_array("deaths_cumulative",
    type = "int",
    dims = list(N, G),
    lower = 0
  )

  # Observation model parameters phi_inv
  phi_D_inv_var <- stan_var("phi_D_inv", lower = 0)
  phi_D_var <- stan_var("phi_D", lower = 0)
  phi_D_inv <- stan_param(phi_D_inv_var, "exponential(5);")
  phi_D <- stan_transform(phi_D_var, "parameters", "inv(phi_D_inv);")

  # All loglik variables
  loglik_vars <- list(delta, deaths_cumulative, phi_D_inv, phi_D)

  # Function bodies
  odefun_body <- "
    vector[4*G] dy_dt; // first G are S, next G are I, etc...
    vector[G] lambda = rep_vector(0.0, G); // force of infection
    for(g in 1:G) {
      for(h in 1:G) {
        lambda[g] += contacts[g,h] * y[G+h]/pop_sizes[h];
      }
    }
    for(g in 1:G) {
      real S_g = y[g];
      real I_g = y[G+g];
      real R_g = y[2*G+g];
      real D_g = y[3*G+g];
      dy_dt[g] = -  beta * lambda[g] * S_g;
      dy_dt[G+g] =  beta * lambda[g] * S_g - (gamma[g] + mu[g]) * I_g;
      dy_dt[2*G+g] = gamma[g]*I_g;
      dy_dt[3*G+g] = mu[g]*I_g;
    }
    return dy_dt;
  "

  loglik_body <- "
    real log_lik = 0.0;
    for(n in 1:N) {
      for(g in 1:G) {
        real D_sol = y_sol[n][3*G + g];
        log_lik += neg_binomial_2_lpmf(deaths_cumulative[n,g] | D_sol + delta,
          phi_D);
      }
    }
    return(log_lik);
  "

  if (prior_only) {
    loglik_body <- ""
    loglik_vars <- list(delta, phi_D_inv, phi_D)
  }

  # Return
  ode_model(
    N = N,
    odefun_vars = odefun_vars,
    odefun_body = odefun_body,
    odefun_init = y0,
    loglik_vars = loglik_vars,
    loglik_body = loglik_body,
    ...
  )
}
