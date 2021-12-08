
# Data
setup_standata_tmdd <- function() {
  t <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0)
  N <- length(t)
  L0 <- 14.8
  data_list <- list(
    N = N,
    t = t,
    L0 = L0,
    D = 3,
    y = rep(1.0, N) # dummy
  )
  return(data_list)
}

# Stan code parts
setup_stancode_tmdd <- function(solver = "rk45") {
  pars <- "
  vector<lower=0>[6] k; // on, off, in, out, eL, eP
  real<lower=0> sigma;
  "
  tpars <- "  
    real R0 = k[3]/k[4];
    vector[3] x0 = to_vector({L0, R0, 0.0});
  "
  prior <- "
  k[1] ~ lognormal(-1, 0.3);  // k_on
  k[2] ~ lognormal(0, 0.3);  // k_off
  k[3] ~ lognormal(0, 0.3);  // k_in
  k[4] ~ lognormal(0, 0.3);  // k_out
  k[5] ~ lognormal(-1, 0.3);  // k_eL
  k[6] ~ lognormal(-3, 0.3);  // k_eP
  sigma ~ lognormal(1, 0.3);
  "
  funs <- "
  // TMDD system right-hand side
  vector TMDD(real t, vector y, vector k) {
    vector[3] dy_dt; // L, R, P
    real L = y[1];
    real R = y[2];
    real P = y[3];
    real rem = k[1]*L*R - k[2]*P;
    dy_dt[1] = - k[5]*L - rem;
    dy_dt[2] = k[3] - k[4]*R - rem;
    dy_dt[3] = rem - k[6]*P;
    return dy_dt;
  }
  "
  data <- "
  int<lower=1> N;            // number of time points
  real t[N];                 // time points
  real L0;                   // initial bolus
  "
  tdata <- "
    real t0 = 0.0;
  "
  obsdata <- "  vector[N] y;"
  if (solver == "rk45") {
    odesolve <- "  vector[3] x[N] = ode_rk45_tol(TMDD, x0, t0, t, rel_tol, abs_tol, max_num_steps, k);"
  } else if (solver == "bdf") {
    odesolve <- "  vector[3] x[N] = ode_bdf_tol(TMDD, x0, t0, t, rel_tol, abs_tol, max_num_steps, k);"
  } else {
    stop("Invalid 'solver' argument!")
  }

  likelihood <- "
  for(n in 1:N) {
    y[n] ~ normal(x[n][3], sigma);
  }
  "
  genquant <- "
  vector[N] y_gen;
  vector[N] L_gen;
  real log_lik = 0.0;
  for(n in 1:N) {
    L_gen[n] = x[n][1];
    y_gen[n] = normal_rng(x[n][1], sigma);
    log_lik += normal_lpdf(y[n] | x[n][1], sigma);
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

# Plotting
plot_tmdd <- function(fit, data) {
  y_gen_rvar <- posterior::as_draws_rvars(fit$draws("y_gen"))$y_gen
  df <- create_ribbon_plot_df(y_gen_rvar)
  df$t <- data$t
  cs <- bayesplot::color_scheme_get()
  fill_alpha <- 1.0
  plt_y <- ggplot(df, aes(
    x = t, y = median, ymin = lower1, ymax = upper1
  )) +
    geom_ribbon(alpha = fill_alpha, fill = cs$light_highlight) +
    geom_ribbon(
      alpha = fill_alpha, fill = cs$mid,
      aes(ymin = lower2, ymax = upper2)
    ) +
    geom_line(col = cs$mid_highlight, lwd = 1) +
    ylab("L")
  if (!is.null(data$y)) {
    df <- data.frame(data$t, data$y)
    colnames(df) <- c("t", "y")
    plt_y <- plt_y + geom_point(
      data = df, aes(x = t, y = y),
      inherit.aes = FALSE, pch = 16
    )
  }
  return(plt_y)
}

# Adding simulated data
add_simulated_data_tmdd <- function(fit, data) {
  y_gen_arr <- posterior::merge_chains(fit$draws("y_gen"))
  S <- dim(y_gen_arr)[1]
  idx <- sample.int(n = S, size = 1)
  y_gen <- subset_draws(y_gen_arr, draw = idx)
  y_gen <- mean(as_draws_rvars(y_gen)$y_gen)
  data$y <- y_gen
  list(data = data, idx = idx)
}
