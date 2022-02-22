# Setup for the Lotka-Volterra experiments

# Following https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html
load_data_lynxhare <- function() {
  df <- odemodeling::lynxhare
  N <- length(df$year) - 1
  ts <- 1:N
  y_init <- c(df$hare[1], df$lynx[1])
  y <- as.matrix(df[2:(N + 1), 2:3])
  y <- cbind(y[, 2], y[, 1]) # hare, lynx order
  colnames(y) <- c("hare", "lynx")

  # Return
  list(
    t0 = df$year[1],
    t = df$year[2:(N + 1)],
    y_obs_init = y_init,
    y_obs = y
  )
}

# Create Lotka-Volterra model
lv_model <- function(prior_only = FALSE, ...) {

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

# Create model and load data
model <- lv_model()
dat <- load_data_lynxhare()
add_data <- list(
  y_obs_init = dat$y_obs_init,
  y_obs = dat$y_obs,
  D = 2
)

# Define initial params for sampling
init <- list(
  alpha = 1,
  beta = 0.1,
  gamma = 1,
  delta = 0.1,
  y0 = dat$y_obs_init,
  sigma = c(1, 1)
)
init <- rep(list(init), CHAINS) # same for all chains
