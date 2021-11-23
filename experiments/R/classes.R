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
  hmc_initial_step_size = 1,
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
    cat("  Pars: {", paste(self$param_names, collapse = ", "), "}\n", sep = "")
    cat("  Solver: ", self$solver, "\n", sep = "")
    cat("  Init is NULL: ", is.null(self$init), "\n", sep = "")
    cat("  HMC initial step size: ", self$hmc_initial_step_size, "\n", sep = "")
    invisible(self)
  },
  sample_prior = function(...) {
    stan_opts <- self$stan_opts
    self$stanmodels$prior$sample(
      data = self$data,
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
  time_posterior_sampling = function(res_dir, idx,
                                     tols, max_num_steps, chains = 4, ...) {
    L <- length(tols)
    WT <- matrix(0.0, L, chains)
    ST <- matrix(0.0, L, chains)
    TT <- matrix(0.0, L, chains)
    GT <- rep(0.0, L)
    j <- 0
    for (tol_j in tols) {
      cat(" (", j, ") Timing posterior sampling with tol = ",
        tol_j, "...\n",
        sep = ""
      )
      j <- j + 1
      fn <- file.path(res_dir, paste0("fit_t_", idx, "_", j, ".rds"))
      sargs <- list(
        abs_tol = tol_j,
        rel_tol = tol_j,
        max_num_steps = max_num_steps
      )
      post_fit <- self$sample_posterior(sargs,
        chains = chains,
        step_size = self$hmc_initial_step_size
      )
      cat("Saving fit to ", fn, "\n", sep = "")
      post_fit$save_object(fn)
      t <- post_fit$time()$chains$total
      gt <- post_fit$time()$total
      GT[j] <- gt
      cat("Grand total time = ", gt, "seconds. \n")
      WT[j, ] <- post_fit$time()$chains$warmup
      ST[j, ] <- post_fit$time()$chains$sampling
      TT[j, ] <- t
    }
    times <- list(warmup = WT, sampling = ST, total = TT, grand_total = GT)
    return(times)
  },
  plot = function(fit) {
    eval(call(paste0("plot_", self$name), fit, self$data))
  },
  add_simulated_data = function(sim) {
    sd <- eval(call(
      paste0("add_simulated_data_", self$name), sim, self$data
    ))
    self$data <- sd$data
    self$sim_par_idx <- sd$idx
  },
  set_init = function(param_draws) {
    self$init <- 0
  },
  set_hmc_initial_step_size = function(step_size) {
    self$hmc_initial_step_size <- step_size
  }
))
