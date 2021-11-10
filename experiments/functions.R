# CLASSES -----------------------------------------------------------------
OdeExperimentSetup <- R6Class("OdeExperimentSetup", list(
  name = NULL,
  stanmodels = NULL,
  solver_args_gen = NULL,
  solver = NULL,
  data = NULL,
  stan_opts = NULL,
  param_names = NULL,
  init = NULL,
  initialize = function(name, solver, solver_args_gen, stan_opts, param_names) {
    self$data <- eval(call(paste0("setup_standata_", name)))
    self$name <- name
    code_parts <- eval(call(paste0("setup_stancode_", name), solver))
    self$stanmodels <- create_cmdstan_models(code_parts)
    self$solver_args_gen <- solver_args_gen
    self$solver <- solver
    self$stan_opts <- stan_opts
    self$param_names <- param_names
  },
  print = function(...) {
    cat("OdeExperimentSetup: \n")
    cat("  Name: ", self$name, "\n", sep = "")
    cat("  Params: ", paste(self$param_names, collapse = ", "), "\n", sep = "")
    cat("  Solver: ", self$solver, "\n", sep = "")
    cat("  Init is NULL: ", is.null(self$init), "\n", sep = "")
    invisible(self)
  },
  sample_prior = function(...) {
    stan_opts <- self$stan_opts
    self$stanmodels$prior$sample(
      sig_figs = stan_opts$sig_figs,
      seed = stan_opts$seed,
      ...
    )
  },
  sample_posterior = function(solver_args, ...) {
    stan_opts <- self$stan_opts
    data <- self$data
    init <- self$init
    sample_posterior(self$stanmodels$posterior, data, solver_args, stan_opts,
      init = init, ...
    )
  },
  time_posterior_sampling = function(tols, max_num_steps = 1e6, chains = 4, ...) {
    L <- length(tols)
    WT <- matrix(0.0, L, chains)
    ST <- matrix(0.0, L, chains)
    TT <- matrix(0.0, L, chains)
    j <- 0
    cat("Timing...\n")
    for (tol_j in tols) {
      cat(" * tol = ", tol_j)
      j <- j + 1
      sargs <- list(
        abs_tol = tol_j,
        rel_tol = tol_j,
        max_num_steps = max_num_steps
      )
      post_fit <- setup$sample_posterior(sargs,
        chains = chains,
        refresh = 0
      )
      t <- post_fit$time()$chains$total
      cat(", t_mean = ", mean(t), ", t_std = ", stats::sd(t), " \n")
      WT[j, ] <- post_fit$time()$chains$warmup
      ST[j, ] <- post_fit$time()$chains$sampling
      TT[j, ] <- t
    }
    times <- list(warmup = WT, sampling = ST, total = TT)
    return(times)
  },
  plot = function(fit) {
    eval(call(paste0("plot_", self$name), fit, self$data))
  },
  add_simulated_data = function(sim) {
    self$data <- eval(call(
      paste0("add_simulated_data_", self$name), sim,
      self$data
    ))
  },
  set_init = function(param_draws) {
    S <- dim(param_draws)[1]
    idx <- sample.int(n = S, size = 1)
    par_init <- param_draws[idx, , ]
    nc <- posterior::nchains(param_draws)
    pin <- list()
    for (idx_c in 1:nc) {
      a <- as.list(par_init[, idx_c, ])
      names(a) <- self$param_names
      pin[[idx_c]] <- a
    }
    self$init <- pin
  }
))


# WORKFLOW ----------------------------------------------------------------

# Run the workflow
run_workflow <- function(setup, tol_init = 1e-4, tol_reduce_factor = 2,
                         max_num_steps = 1e6) {

  # Sample posterior using tol_init
  sargs_sample <- list(
    abs_tol = tol_init,
    rel_tol = tol_init,
    max_num_steps = max_num_steps
  )
  cat("Sampling posterior with tol_init (", tol_init, ")...\n", sep = "")
  post_fit <- setup$sample_posterior(
    solver_args = sargs_sample,
    refresh = 0
  )
  post_draws <- post_fit$draws(setup$param_names)
  cat("Done.\n")

  # Create a plot using posterior draws
  post_sim <- simulate(setup, post_draws, sargs_sample)
  plot_sim_plot <- setup$plot(post_sim)

  # Tune the reference method so that it is reliable at post_draws
  tuning <- tune_solver_tols(
    setup, post_sim, post_draws, sargs_sample, tol_reduce_factor
  )
  tuning_plot <- plot_tuning(tuning)

  # Compute importance weights and Pareto-$k$
  is <- use_psis(post_sim, tuning$last_sim)

  # Importance resampling
  post_draws_resampled <- posterior::resample_draws(
    post_draws,
    weights = is$log_weights
  )

  # Timing
  workflow_time <- post_fit$time()$total + post_sim$time()$total + tuning$total_time
  cat("Total workflow time was", workflow_time, "seconds.\n", sep = " ")

  # Return
  list(
    post_fit = post_fit,
    post_sim = post_sim,
    post_sim_plot = post_sim_plot,
    tuning = tuning,
    tuning_plot = tuning_plot,
    is = is,
    post_draws = post_draws,
    post_draws_resampled = post_draws_resampled,
    workflow_time = workflow_time
  )
}

# UTILS -------------------------------------------------------------------

# Validate solver arguments
check_sa <- function(solver_args) {
  MAX_INT <- 2^31 - 1
  required <- c("rel_tol", "abs_tol", "max_num_steps")
  checkmate::assert_names(names(solver_args), permutation.of = required)
  checkmate::assertNumeric(solver_args$rel_tol, lower = 0)
  checkmate::assertNumeric(solver_args$abs_tol, lower = 0)
  checkmate::assertIntegerish(solver_args$max_num_steps, lower = 1, upper = MAX_INT)
  TRUE
}

# Print output if running Stan model failed
print_output_if_failed <- function(stan_out) {
  codes <- stan_out$return_codes()
  idx_failed <- which(codes > 0)
  for (idx in idx_failed) {
    cat("Chain ", idx, ", failed, printing its output:\n", sep = "")
    print(stan_out$output(idx))
  }
  if (length(idx_failed == 0)) {
    cat("All chains were successful.\n")
  }
}


# CREATING STAN MODELS ------------------------------------------------------

# Create a block of Stan code
stan_block <- function(block_name, block_code) {
  block_name <- trimws(tolower(block_name))
  block_code <- trimws(block_code, whitespace = "[\n]")
  paste0(block_name, " {\n", block_code, "\n}\n")
}

# Create many blocks of Stan code
stan_code <- function(blocks, codes) {
  model_code <- ""
  J <- length(blocks)
  for (j in seq_len(J)) {
    model_code <- paste(model_code, stan_block(blocks[j], codes[j]), sep = "\n")
  }
  return(model_code)
}

# Create full Stan code for prior sampling model
prior_model_code <- function(pars, tpars, prior) {
  cat("* Creating CmdStanModel for prior sampling...\n")
  blocks <- c("parameters", "transformed parameters", "model")
  codes <- c(pars, tpars, prior)
  stan_code(blocks, codes)
}

# Create full Stan code for simulator model
simulator_model_code <- function(funs, data, tdata, pars, tpars, ode, gq, od) {
  cat("* Creating CmdStanModel for simulating...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "generated quantities"
  )
  gq <- paste(ode, gq, sep = "\n")
  data <- paste(data, od, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, gq)
  stan_code(blocks, codes)
}

# Create full Stan code for posterior sampling model
posterior_model_code <- function(funs, data, tdata, obsdata,
                                 pars, tpars, prior, ode, lik) {
  cat("* Creating CmdStanModel for posterior sampling...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "model"
  )
  post <- paste(prior, ode, lik, sep = "\n")
  data <- paste(data, obsdata, sep = "\n")
  codes <- c(funs, data, tdata, pars, tpars, post)
  stan_code(blocks, codes)
}

# Create all CmdStanModels
create_cmdstan_models <- function(code) {
  pars <- code$pars
  tpars <- code$tpars
  prior <- code$prior
  funs <- code$functions
  ode <- code$odesolve
  obsdata <- code$obsdata
  lik <- code$likelihood
  data <- code$data
  tdata <- code$tdata
  gq <- code$genquant
  codes <- list(
    prior = prior_model_code(pars, tpars, prior),
    simulator = simulator_model_code(
      funs, data, tdata, pars, tpars, ode, gq, obsdata
    ),
    posterior = posterior_model_code(
      funs, data, tdata, obsdata,
      pars, tpars, prior, ode, lik
    )
  )
  j <- 0
  models <- list()
  for (code in codes) {
    j <- j + 1
    models[[j]] <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  }
  names(models) <- names(codes)
  return(models)
}


# RUNNING CMDSTAN -----------------------------------------------------

# Function for simulating ODE solutions and data given parameter(draws)s
simulate <- function(setup, params, solver_args) {
  stopifnot(is(setup, "OdeExperimentSetup"))
  stopifnot(is(params, "draws"))
  stopifnot(is(solver_args, "list"))
  check_sa(solver_args)
  data <- setup$data
  stan_opts <- setup$stan_opts
  model <- setup$stanmodels$simulator
  capture.output({
    out <- model$generate_quantities(
      data = c(data, solver_args),
      fitted_params = params,
      seed = stan_opts$seed,
      sig_figs = stan_opts$sig_figs
    )
  })
  print_output_if_failed(out)
  return(out)
}

# Using simulate with different tolerances
simulate_many <- function(setup, params, atol, rtol) {
  stopifnot(is(params, "draws"))
  data <- setup$data
  max_num_steps <- setup$solver_args_gen$max_num_steps
  J1 <- length(atol)
  J2 <- length(rtol)
  S <- posterior::niterations(params) * posterior::nchains(params)
  TIME <- array(0, dim = c(J1, J2))
  XSIM <- array(0, dim = c(J1, J2, S, data$N * data$D))
  LL <- array(0, dim = c(J1, J2, S))
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      solver_args <- list(
        abs_tol = atol[j1],
        rel_tol = rtol[j2],
        max_num_steps = max_num_steps
      )
      tryCatch(
        expr = {
          sim <- simulate(setup, params, solver_args)
          XSIM[j1, j2, , ] <- posterior::merge_chains(
            sim$draws("x")
          )[, 1, , drop = TRUE]
          TIME[j1, j2] <- sim$time()$total
          LL[j1, j2, ] <- posterior::merge_chains(
            sim$draws("log_lik")
          )[, 1, 1, drop = TRUE]
        },
        error = function(e) {
          message(
            "Caught an error in simulate! Likely max_num_steps was too",
            " low or negative solution was obtained so likelihood could",
            " not be computed!"
          )
          stop(e)
        }
      )
    }
  }
  return(list(times = TIME, sims = XSIM, log_liks = LL))
}

# Function for posterior sampling
sample_posterior <- function(model, data, solver_args, stan_opts, ...) {
  stopifnot(is(model, "CmdStanModel"))
  stopifnot(is(data, "list"))
  stopifnot(is(solver_args, "list"))
  check_sa(solver_args)
  fit <- model$sample(
    data = c(data, solver_args),
    sig_figs = stan_opts$sig_figs,
    seed = stan_opts$seed,
    ...
  )
  print_output_if_failed(fit)
  return(fit)
}


# COMPUTING ERRORS --------------------------------------------------------

# Compute error in x compared to x_ref
compute_sol_error <- function(x, x_ref, fun) {
  abs_errors <- abs(as.vector(x_ref) - as.vector(x))
  eval(call(fun, abs_errors))
}

# Compute error to most accurate solution
compute_sol_errors <- function(XSIM, fun = "max") {
  J1 <- dim(XSIM)[1]
  J2 <- dim(XSIM)[2]
  ERR <- array(0, dim = c(J1, J2))
  x_ref <- XSIM[1, 1, , ]
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      x <- XSIM[j1, j2, , ]
      ERR[j1, j2] <- compute_sol_error(x, x_ref, fun)
    }
  }
  return(ERR)
}

# Compute error in likelihood, compared to most accurate solution
compute_loglik_errors <- function(LL, fun = "max") {
  J1 <- dim(LL)[1]
  J2 <- dim(LL)[2]
  ERR <- array(0, dim = c(J1, J2))
  loglik_best <- LL[1, 1, ]
  for (j1 in 1:J1) {
    for (j2 in 1:J2) {
      loglik <- LL[j1, j2, ]
      abs_errors <- abs(as.vector(loglik_best) - as.vector(loglik))
      ERR[j1, j2] <- eval(call(fun, abs_errors))
    }
  }
  return(ERR)
}


# COMPUTING PARETO_K --------------------------------------------------------

# Compute log importance weights
log_importance_weights <- function(fit_high, fit_low) {
  stopifnot(is(fit_high, "CmdStanFit"))
  stopifnot(is(fit_low, "CmdStanFit"))
  LL_high <- fit_high$draws("log_lik")[, , 1, drop = TRUE]
  LL_low <- fit_low$draws("log_lik")[, , 1, drop = TRUE]
  return(LL_high - LL_low)
}

# Compute importance weights and pareto k diagnostic
use_psis <- function(fit_high, fit_low) {
  log_ratios <- log_importance_weights(fit_high, fit_low)
  chain_id <- rep(1:ncol(log_ratios), each = nrow(log_ratios))
  x <- as.vector(exp(-log_ratios))
  r_eff <- loo::relative_eff(x, chain_id)
  loo::psis(x, r_eff = r_eff)
}


# TUNING THE SOLVER -------------------------------------------------------

# Extract ODE solutions
get_x_sim <- function(sim) {
  posterior::merge_chains(sim$draws("x"))[, 1, , drop = TRUE]
}

# Run tuning
tune_solver_tols <- function(setup, p_sim, p_params, p_sargs, factor) {
  S <- posterior::ndraws(posterior::merge_chains(p_sim$draws()))
  tol_j <- min(p_sargs$rel_tol, p_sargs$abs_tol)
  data <- setup$data
  stan_opts <- setup$stan_opts
  model <- setup$stanmodels$simulator
  error_to_ref <- Inf
  p_x <- get_x_sim(p_sim)
  p_t <- p_sim$time()$total
  res <- c(1 / tol_j, p_t, 0.0, NA, 1.0)
  idx <- 0
  cat("Tuning tolerances...\n")
  while (idx < 13) {
    idx <- idx + 1
    tol_j <- tol_j / factor
    cat(" * tol = ", tol_j, sep = "")
    sargs <- list(
      abs_tol = tol_j,
      rel_tol = tol_j,
      max_num_steps = 10^9
    )
    sim <- simulate(setup, p_params, sargs)
    stopifnot(is(sim, "CmdStanFit"))
    t_j <- sim$time()$total
    cat(", time = ", t_j, " s", sep = "")
    x <- get_x_sim(sim)
    err_j <- compute_sol_error(x, p_x, "max")
    cat(", mae = ", err_j, sep = "")
    is <- use_psis(sim, p_sim)
    k_j <- is$diagnostics$pareto_k
    n_j <- is$diagnostics$n_eff / S
    cat(", k_hat = ", k_j, ", n_eff = ", n_j, "\n", sep = "")
    res_j <- c(1 / tol_j, t_j, err_j, k_j, n_j)
    res <- rbind(res, res_j)
  }
  colnames(res) <- c("inv_tol", "time", "mae", "k_hat", "p_eff")
  res <- data.frame(res)
  rownames(res) <- NULL

  # Return
  list(
    metrics = res,
    last_sim = sim,
    total_time = sum(res$time),
    max_khat = max(res$k_hat, na.rm = TRUE)
  )
}

# Plot tuning results
plot_tuning <- function(tuning, ...) {
  df <- tuning$metrics
  add_geoms <- function(x) {
    x + geom_line() + geom_point() + scale_x_log10()
  }
  p_A <- add_geoms(ggplot(df, aes(x = inv_tol, y = mae)))
  p_B <- add_geoms(ggplot(df, aes(x = inv_tol, y = k_hat)))
  p_C <- add_geoms(ggplot(df, aes(x = inv_tol, y = p_eff)))
  p_D <- add_geoms(ggplot(df, aes(x = inv_tol, y = time)))
  plt <- ggpubr::ggarrange(p_A, p_B, p_C, p_D, labels = "auto", ...)
  return(plt)
}

# PLOTTING ----------------------------------------------------------------

# Plotting helper
create_ribbon_plot_df <- function(rvar) {
  alpha1 <- 0.1
  alpha2 <- 0.25
  c1 <- 100 * (1 - 2 * alpha1)
  c2 <- 100 * (1 - 2 * alpha2)
  cat("Plotting median and central ", c1, "% and ", c2, "% intervals.\n",
    sep = ""
  )
  lower1 <- as.vector(quantile(rvar, probs = alpha1))
  upper1 <- as.vector(quantile(rvar, probs = 1 - alpha1))
  lower2 <- as.vector(quantile(rvar, probs = alpha2))
  upper2 <- as.vector(quantile(rvar, probs = 1 - alpha2))
  median <- as.vector(quantile(rvar, probs = 0.5))
  df <- data.frame(median, lower1, upper1, lower2, upper2)
  return(df)
}
