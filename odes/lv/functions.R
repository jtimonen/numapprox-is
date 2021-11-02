# A function that runs the inference
run_inference <- function(model, data, ITER, CHAINS, ADAPT_DELTA) {

  # Run sampling
  fit <- sampling(
    object = model,
    data = data,
    iter = ITER,
    chains = CHAINS,
    control = list(adapt_delta = ADAPT_DELTA),
    init = "0",
    save_warmup = FALSE
  )

  # Extract log posterior values (not Jacobian adjusted)
  lh1 <- get_samples(fit, "log_lik_na")
  lh2 <- get_samples(fit, "log_lik_na_REF_")
  pr1 <- get_samples(fit, "log_prior_na")
  pr2 <- get_samples(fit, "log_prior_na_REF_")
  post1 <- as.vector(lh1) + as.vector(pr1)
  post2 <- as.vector(lh2) + as.vector(pr2)

  # PSIS
  out <- psis(post1 - post2)
  pareto_k <- out$diagnostics$pareto_k
  runtimes <- rowSums(get_elapsed_time(fit))

  # Return list
  return(list(pareto_k = pareto_k, runtimes = runtimes, fit = fit))
}


# Create additional Stan data related to some ODE solver methods
add_interpolation_data <- function(data_list, h) {
  t0 <- data_list$t0
  ts <- data_list$ts
  R <- compute_R(t0, ts, h)
  A <- compute_A(t0, ts, h, R)
  data_list$STEP_SIZE <- h
  data_list$INTERP_R <- R
  data_list$INTERP_A <- pmax(A, 0) # rounding errors can make A slightly negative though should be 0
  return(data_list)
}

# Compute integers r_1, ..., r_n
compute_R <- function(t0, ts, h) {
  n <- length(ts)
  R <- rep(0, n)
  for (i in 1:n) {
    r <- 0
    while (t0 + r * h < ts[i]) {
      r <- r + 1
    }
    R[i] <- r - 1
  }
  return(R)
}

# Compute multipliers a_1, ..., a_n
compute_A <- function(t0, ts, h, R) {
  n <- length(ts)
  A <- rep(0, n)
  for (i in 1:n) {
    D_i <- ts[i] - (t0 + R[i] * h)
    A[i] <- (h - D_i) / h
  }
  return(A)
}

# Helper function
get_samples <- function(stan_fit, param) {
  samples <- rstan::extract(stan_fit, pars = param)[[param]]
  return(samples)
}
