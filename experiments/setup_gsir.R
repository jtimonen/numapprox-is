# Example contact matrix
create_contact_matrix <- function(G) {
  M <- diag(G) + matrix(runif(G * G), G, G)
  M <- M + t(M)
  return(M / G)
}

# Data
setup_standata_gsir <- function() {
  t <- as.numeric(seq(0.2, 15, by = 0.2))
  N <- length(t)
  G <- 4 # number of groups
  data_list <- list(
    N = N,
    t = t,
    I0 = c(5, 0, 0, 0),
    pop_sizes = c(2000, 4000, 3000, 1000),
    contacts = create_contact_matrix(G),
    D = 2 * G,
    G = G,
    I_data = matrix(1.0, N, G) # dummy
  )
  return(data_list)
}

# Stan code parts
setup_stancode_gsir <- function(solver = "rk45") {
  pars <- "
    real<lower=0> beta;
    vector<lower=0>[G] gamma;
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
  vector SIR(real t, vector y, real beta, vector gamma, data matrix contacts,
      data vector pop_sizes) {
    int G = num_elements(pop_sizes);
    vector[2*G] dy_dt; // first G are susceptible, next G are infected
    vector[G] infection_rates;
    vector[G] recovery_rates;
    vector[G] lambda = rep_vector(0.0, G);
    for(g in 1:G){
      for(h in 1:G) {
        lambda[g] += contacts[g,h] * y[G+h]/pop_sizes[h];
      }
    }
    for(g in 1:G){
      dy_dt[g] = -  beta * lambda[g] * y[g];
      dy_dt[G+g] =  beta * lambda[g] * y[g] - gamma[g] * y[G+g];
    }
    return dy_dt;
  }
  "
  data <- "
    int<lower=1> N;                 // number of time points
    real t[N];                      // time points
    int<lower=1> G;                 // number of groups
    vector[G] pop_sizes;            // population sizes in each group
    vector[G] I0;                   // initial number of infected in each group
    matrix[G, G] contacts;          // contact matrix
    real<lower=0> rel_tol;          // ODE solver relative tolerance
    real<lower=0> abs_tol;          // ODE solver absolute tolerance
    int<lower=0> max_num_steps;     // ODE solver maximum number of steps
  "
  tdata <- "
    real t0 = 0.0;
    vector[2*G] x0;
    for(g in 1:G){
      x0[g] = pop_sizes[g] - I0[g]; // S
    }
    for(g in 1:G){
      x0[G + g] = I0[g]; // I
    }
  "
  obsdata <- "  int<lower=0> I_data[N, G];"
  if (solver == "rk45") {
    odesolve <- "  vector[2*G] x[N] = ode_rk45_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, beta, gamma, contacts, pop_sizes);"
  } else if (solver == "bdf") {
    odesolve <- "  vector[2*G] x[N] = ode_bdf_tol(SIR, x0, t0, t, rel_tol, abs_tol, max_num_steps, beta, gamma, contacts, pop_sizes);"
  } else {
    stop("Invalid 'solver' argument!")
  }

  likelihood <- "
    for(n in 1:N) {
      for(g in 1:G) {
        I_data[n,g] ~ neg_binomial_2(x[n][G+g] + 10*abs_tol, phi);
      }
    }
  "
  genquant <- "
    int I_gen[N, G];
    real log_lik = 0.0;
    for(n in 1:N) {
      for(g in 1:G) {
        I_gen[n,g] = neg_binomial_2_rng(x[n][G+g] + 10*abs_tol, phi);
        log_lik += neg_binomial_2_lpmf(I_data[n,g] | x[n][G+g] + 10*abs_tol, phi);
      }
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
plot_gsir <- function(fit, data) {
  I_gen_rvar <- posterior::as_draws_rvars(fit$draws("I_gen"))$I_gen
  G <- data$G
  df <- NULL
  for (g in 1:G) {
    df_g <- create_ribbon_plot_df(I_gen_rvar[, g])
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
add_simulated_data_gsir <- function(fit, data) {
  I_gen_arr <- posterior::merge_chains(fit$draws("I_gen"))
  S <- dim(I_gen_arr)[1]
  idx <- sample.int(n = S, size = 1)
  I_gen <- subset_draws(I_gen_arr, draw = idx)
  I_gen <- mean(as_draws_rvars(I_gen)$I_gen)
  data$I_data <- I_gen
  list(data = data, idx = idx)
}
