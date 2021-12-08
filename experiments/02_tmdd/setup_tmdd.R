
# Data
setup_standata_tmdd <- function() {
  t <- as.numeric(seq(1, 15, by = 0.2))
  N <- length(t)
  L0 <- 1
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
  k ~ normal(0.3, 0.5);
  sigma ~ normal(1, 0.1);
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
  real L0;                   // initial number of infected in each group
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
  real log_lik = 0.0;
  for(n in 1:N) {
    y_gen[n] = normal_rng(x[n][3], sigma);
    log_lik += normal_lpdf(y[n] | x[n][3], sigma);
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
  G <- data$G
  df <- NULL
  for (g in 1:G) {
    df_g <- create_ribbon_plot_df(y_gen_rvar[, g])
    df_g$group <- rep(g, nrow(df_g))
    df <- rbind(df, df_g)
  }
  df$group <- as.factor(df$group)

  df$t <- data$t
  cs <- bayesplot::color_scheme_get()
  fill_alpha <- 1.0
  plt_I <- ggplot(df, aes(
    x = t, y = median, ymin = lower1, ymax = upper1,
    group = group
  )) +
    geom_ribbon(alpha = fill_alpha, fill = cs$light_highlight) +
    geom_ribbon(
      alpha = fill_alpha, fill = cs$mid,
      aes(ymin = lower2, ymax = upper2)
    ) +
    geom_line(col = cs$mid_highlight, lwd = 1) +
    ylab("I") +
    facet_wrap(. ~ group)
  if (!is.null(data$I_data)) {
    df <- data.frame(data$t, data$I_data)
    colnames(df) <- c("t", 1:G)
    df_long <- pivot_longer(df, cols = as.character(1:G))
    colnames(df_long) <- c("t", "group", "I_data")
    df_long$group <- as.factor(df_long$group)
    plt_I <- plt_I + geom_point(
      data = df_long, aes(x = t, y = I_data, group = group),
      inherit.aes = FALSE, pch = 16
    )
  }
  return(plt_I)
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
