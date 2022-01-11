# Validate fit
validate_fit <- function(setup, sampled, tols, max_num_steps, ...) {
  post_fit <- sampled$fit
  post_draws <- post_fit$draws(setup$param_names)
  sargs <- list(
    abs_tol = sampled$tol,
    rel_tol = sampled$tol,
    max_num_steps = max_num_steps
  )
  post_sim <- simulate(setup, post_draws, sargs)
  
  # Simulations
  sims_reference <- simulate_many(setup, post_draws, tols, max_num_steps)
  times <- sapply(sims_reference, function(x) {
    x$time()$total
  })
  
  # Tune the reference method so that it is reliable at post_draws
  cat("\nValidating tolerances...\n")
  tuning <- validate_tols(
    setup, post_sim, sims_reference
  )
  # tuning_plot <- plot_tuning(tuning)
  list(
    time = post_sim$time()$total,
    sims_reference = sims_reference,
    times_reference = times,
    tuning = tuning
  )
}
