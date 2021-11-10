# Data
setup_standata_sir <- function() {
  t <- as.numeric(seq(0.2, 15, by = 0.2))
  N <- length(t)
  data_list <- list(
    N = N,
    t = t,
    I0 = 5,
    pop_size = 1000,
    D = 2,
    I_data = rep(1.0, N) # dummy
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
    gamma ~ normal(0.3, 0.3);
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
    int I_gen[N];
    real log_lik = 0.0;
    for(n in 1:N) {
      I_gen[n] = neg_binomial_2_rng(x[n][2] + 10*abs_tol, phi);
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

# Plotting
plot_sir <- function(fit, data) {
  I_gen_rvar <- posterior::as_draws_rvars(fit$draws("I_gen"))$I_gen
  df <- create_ribbon_plot_df(I_gen_rvar)
  df$t <- data$t
  cs <- bayesplot::color_scheme_get()
  fill_alpha <- 1.0
  plt_I <- ggplot(df, aes(x = t, y = median, ymin = lower1, ymax = upper1)) +
    geom_ribbon(alpha = fill_alpha, fill = cs$light_highlight) +
    geom_ribbon(
      alpha = fill_alpha, fill = cs$mid,
      aes(ymin = lower2, ymax = upper2)
    ) +
    geom_line(col = cs$mid_highlight, lwd = 1) +
    ylab("I")
  if (!is.null(data$I_data)) {
    df <- data.frame(data$t, data$I_data)
    colnames(df) <- c("t", "I_data")
    plt_I <- plt_I + geom_point(
      data = df, aes(x = t, y = I_data),
      inherit.aes = FALSE, pch = 16
    )
  }
  return(plt_I)
}

# Adding simulated data
add_simulated_data_sir <- function(fit, data) {
  I_gen_arr <- posterior::merge_chains(fit$draws("I_gen"))
  S <- dim(I_gen_arr)[1]
  idx <- sample.int(n = S, size = 1)
  data$I_data <- as.vector(I_gen_arr[idx, 1, ])
  list(data = data, idx = idx)
}
