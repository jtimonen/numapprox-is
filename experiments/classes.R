# R6 class
OdeExperimentSetup <- R6Class("OdeExperimentSetup", list(
  name = NULL,
  stanmodels = NULL,
  solver_args_gen = NULL,
  solver = NULL,
  data = NULL,
  stan_opts = NULL,
  param_names = NULL,
  init = NULL,
  sim_par_idx = NULL,
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
    cat("  Params: {", paste(self$param_names, collapse = ", "), "}\n", sep = "")
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
      cat(" * Timing posterior sampling with tol = ", tol_j, "...\n", sep = "")
      j <- j + 1
      sargs <- list(
        abs_tol = tol_j,
        rel_tol = tol_j,
        max_num_steps = max_num_steps
      )
      post_fit <- setup$sample_posterior(sargs,
        chains = chains,
        refresh = 0,
        show_messages = FALSE
      )
      t <- post_fit$time()$chains$total
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
    sd <- eval(call(
      paste0("add_simulated_data_", self$name), sim,
      self$data
    ))
    self$data <- sd$data
    self$sim_par_idx <- sd$idx
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
