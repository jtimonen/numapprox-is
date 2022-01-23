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
