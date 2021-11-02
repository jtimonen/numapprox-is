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
simulator_model_code <- function(funs, data, tdata, pars, tpars, ode, gq) {
  cat("* Creating CmdStanModel for simulating...\n")
  blocks <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters", "generated quantities"
  )
  gq <- paste(ode, gq, sep = "\n")
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
      funs, data, tdata, pars, tpars, ode, gq
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

# Function for simulating ODE solutions and data given parameter( draws)s
simulate <- function(model, params, data, solver_args, stan_opts) {
  stopifnot(is(model, "CmdStanModel"))
  stopifnot(is(params, "draws"))
  stopifnot(is(data, "list"))
  stopifnot(is(solver_args, "list"))
  check_sa(solver_args)
  out <- model$generate_quantities(
    data = c(data, solver_args),
    fitted_params = params,
    seed = stan_opts$seed,
    sig_figs = stan_opts$sig_figs
  )
  print_output_if_failed(out)
  return(out)
}

# Using simulate with different tolerances
simulate_many <- function(model, params, data, stan_opts,
                          atol, rtol, max_num_steps) {
  stopifnot(is(params, "draws"))
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
          sim <- simulate(model, params, data, solver_args, stan_opts)
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

tune_solver <- function(tols, model, params, data, stan_opts,
                        max_num_steps, err_tol) {
  error_to_ref <- Inf
  x_ref <- NULL
  TIMES <- c()
  ERR <- c()
  idx <- 0
  for (TOL in tols) {
    idx <- idx + 1
    cat("* simulating with abs_tol = rel_tol = ", TOL, "\n", sep = "")
    sargs <- list(
      abs_tol = TOL,
      rel_tol = TOL,
      max_num_steps = max_num_steps
    )
    sim <- simulate(model, params, data, sargs, stan_opts)
    x <- posterior::merge_chains(sim$draws("x"))[, 1, , drop = TRUE]
    TIMES <- c(TIMES, sim$time()$total)
    if (is.null(x_ref)) {
      x_ref <- x
    }
    err <- compute_sol_error(x, x_ref, "max")
    ERR <- c(ERR, err)
    if (idx > 1 && err < err_tol) {
      return(list(
        tols = tols[1:idx], errors = ERR, times = TIMES, last_sim
        = sim, last_tol = tols[idx]
      ))
    }
    x_ref <- x
  }
  warning("err_tol was not reached")
  list(
    tols = tols[1:idx], errors = ERR, times = TIMES, last_sim = sim,
    last_tol = tols[idx]
  )
}


# PLOTTING ----------------------------------------------------------------

# Runtimes plot
plot_sim_times <- function(abs_tol, rel_tol, TIME) {
  par(mfrow = c(1, 2))
  plot(log10(rel_tol), diag(TIME),
    xlab = "log10(tol)", ylab = "time (s)",
    type = "o", pch = 16
  )
  grid()
  image(log10(abs_tol), log10(rel_tol), TIME, main = "time (s)")
}


# Errors plot
plot_sim_errors <- function(abs_tol, rel_tol, ERR, log = TRUE) {
  main <- deparse(substitute(ERR))
  if (log) {
    ERR <- log(ERR)
    main <- paste0("log(", main, ")")
  }
  par(mfrow = c(1, 2))
  plot(log10(rel_tol), diag(ERR),
    xlab = "log10(tol)", ylab = main,
    type = "o", pch = 16
  )
  grid()
  image(log10(abs_tol), log10(rel_tol), ERR, main = main)
}
