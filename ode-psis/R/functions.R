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


# Plotting metrics --------------------------------------------------------

# Plotting helper function
plot_metric.create_df <- function(confs, values, conf_name, metric_name) {
  checkmate::assert_numeric(confs, len = length(values))
  df <- data.frame(confs, values)
  colnames(df) <- c(conf_name, metric_name)
  return(df)
}

# Plotting helper function
plot_metric.create_plot <- function(confs, values, conf_name, metric_name,
                                    log10) {
  df <- plot_metric.create_df(confs, values, conf_name, metric_name)
  plt <- ggplot(df, aes_string(x = conf_name, y = metric_name))
  plt <- plt + geom_line() + geom_point() + theme_bw()
  if (log10) {
    plt <- plt + scale_x_log10(breaks = confs)
  } else {
    plt <- plt + scale_x_continuous(breaks = confs)
  }
  plt + theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank()
  )
}

# Plot a metric
plot_metric <- function(metrics, name, tols = NULL, num_steps = NULL) {
  values <- metrics[, name]
  checkmate::assert_string(name)
  if (!is.null(tols)) {
    if (!is.null(num_steps)) {
      stop("either tols or num_steps must be NULL")
    }
    confs <- 1 / tols
    conf_name <- "inv_tol"
    log10 <- TRUE
  } else {
    if (is.null(num_steps)) {
      stop("either tols or num_steps must not be NULL")
    }
    confs <- num_steps
    conf_name <- "num_steps"
    log10 <- FALSE
  }
  plot_metric.create_plot(confs, values, conf_name, name, log10)
}

# Edit x label
add_inv_tol_xlab <- function(plt) {
  plt + xlab(expression(tol^{
    -1
  }))
}

# Plot pareto_k metric
plot_pareto_k <- function(metrics, tols = NULL, num_steps = NULL) {
  plt <- plot_metric(metrics, "pareto_k", tols, num_steps)
  plt <- plt + geom_hline(yintercept = 0.5, lty = 2, color = "firebrick3")
  plt <- plt + geom_hline(yintercept = 0.7, lty = 2, color = "steelblue")
  plt <- plt + ylab("Pareto-k") 
  add_inv_tol_xlab(plt)
}
