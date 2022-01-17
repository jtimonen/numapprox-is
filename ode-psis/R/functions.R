library(ggplot2)

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

# Plot all metrics
plot_metrics <- function(rel, tols = NULL, num_steps = NULL, flat = FALSE,
                         arrange = TRUE) {
  library(ggpubr)
  if (flat) {
    nrow <- 1
    ncol <- 5
  } else {
    nrow <- 2
    ncol <- 3
  }
  pA <- plot_pareto_k(rel, tols = tols, num_steps = num_steps)
  pB <- plot_r_eff(rel, tols = tols, num_steps = num_steps)
  pC <- plot_metric(rel$times, "time", tols = tols, num_steps = num_steps) +
    ylab("Time (s)")
  pD <- plot_mad(rel, tols = tols, num_steps = num_steps, loglik = FALSE)
  pE <- plot_mad(rel, tols = tols, num_steps = num_steps, loglik = TRUE)
  plots <- list(pA, pB, pC, pD, pE)
  if (arrange) {
    plots <- ggpubr::ggarrange(
      plotlist = plots, labels = "auto",
      nrow = nrow, ncol = ncol
    )
  }
  return(plots)
}
