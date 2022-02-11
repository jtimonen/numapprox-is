library(ggplot2)
library(ggpubr)

# Get relative efficencies
get_diags_df <- function(fits) {
  x <- sapply(fits$files, get_diags)
  colnames(x) <- NULL
  x <- data.frame(t(x))
  colnames(x) <- c("reff", "max_rhat")
  x
}

# Get sampling diagnostics
get_diags <- function(file) {
  fit <- readRDS(file)
  a <- exp(fit$loglik())
  reff <- as.numeric(loo::relative_eff(a))
  max_rhat <- max(fit$summary()[, "rhat"], na.rm = TRUE)
  c(reff, max_rhat)
}

# Load a fit from file and compile
load_fit <- function(file) {
  fit <- readRDS(file)
  ensure_compiled(fit)
}

# Compile model in fit if not already compiled
ensure_compiled <- function(fit) {
  tryCatch(
    {
      fit$model$assert_stanfile_exists()
    },
    error = function(e) {
      fit$model$reinit()
    }
  )
  fit$model$assert_stanfile_exists()
  return(fit)
}

# Plot all metrics
plot_metrics <- function(rel, tols = NULL, num_steps = NULL, flat = FALSE,
                         arrange = TRUE) {
  if (flat) {
    pdims <- c(1, 5)
  } else {
    pdims <- c(2, 3)
  }
  pA <- plot_pareto_k(rel, tols = tols, num_steps = num_steps)
  pareto_k <- rel$metrics[, "pareto_k"]
  ymin <- min(pareto_k)
  ymax <- max(pareto_k)
  if (is.finite(ymax)) {
    pA <- pA + ylim(c(ymin, ymax))
  }
  pB <- plot_r_eff(rel, tols = tols, num_steps = num_steps)
  pC <- plot_metric(rel$times, "time", tols = tols, num_steps = num_steps) +
    ylab("Time (s)")
  pD <- plot_mad(rel, tols = tols, num_steps = num_steps, loglik = FALSE)
  pE <- plot_mad(rel, tols = tols, num_steps = num_steps, loglik = TRUE)
  plots <- list(pA, pB, pC, pD, pE)
  if (arrange) {
    plots <- ggpubr::ggarrange(
      plotlist = plots, labels = "auto",
      nrow = pdims[1], ncol = pdims[2]
    )
  }
  return(plots)
}

# Get tolerances vector
get_tol_vec <- function(solvers) {
  get_tol <- function(x) {
    if (x$abs_tol != x$rel_tol) {
      stop("abs_tol and rel_tol are not equal")
    }
    x$abs_tol
  }
  sapply(solvers, get_tol)
}

# Get num_steps vector
get_num_steps_vec <- function(solvers) {
  get_ns <- function(x) x$num_steps
  sapply(solvers, get_ns)
}

# Plot time comparison (different tolerances)
plot_time_comparison_tol <- function(fits, reliab, idx_ok, ylog = FALSE) {
  df <- create_time_comparison_df(fits, reliab, idx_ok, ylog)
  map <- aes(x = inv_tol, y = time, group = procedure, color = procedure)
  plt <- ggplot(df, map)
  plt <- add_plot_geoms(plt, df$inv_tol, TRUE)
  if (ylog) {
    plt <- plt + ylab("log(time)")
  }
  return(plt)
}

# Plot time comparison (different number of steps)
plot_time_comparison_ns <- function(fits, reliab, idx_ok, ylog = FALSE) {
  df <- create_time_comparison_df(fits, reliab, idx_ok, ylog)
  plt <- ggplot(df, aes(x = num_steps, y = time, group = procedure, color = procedure))
  plt <- add_plot_geoms(plt, df$num_steps, FALSE)
  if (ylog) {
    plt <- plt + ylab("log(time)")
  }
  return(plt)
}

# Create data frame for plotting times
create_time_comparison_df <- function(fits, reliab, idx_ok, ylog) {
  has_tol <- is(fits$solvers[[1]], "AdaptiveOdeSolver")
  if (has_tol) {
    x <- get_tol_vec(fits$solvers)
    x_rel <- get_tol_vec(reliab$solvers)
    x <- c(x, x_rel)
    x <- 1 / x
    x_name <- "inv_tol"
  } else {
    x <- get_num_steps_vec(fits$solvers)
    x_rel <- get_num_steps_vec(reliab$solvers)
    x <- c(x, x_rel)
    x_name <- "num_steps"
  }

  times <- fits$times$grand_total
  t_sample <- times[idx_ok]
  times_rel <- reliab$times + t_sample
  t <- c(times, times_rel)
  if (ylog) {
    t <- log(t)
  }
  l1 <- "high"
  l2 <- "low+psis"
  fac <- c(rep(l1, length(times)), rep(l2, length(times_rel)))
  fac <- factor(x = fac)
  df <- data.frame(fac, t, x)
  colnames(df) <- c("procedure", "time", x_name)
  return(df)
}

# Add geoms and theme to ggplot
add_plot_geoms <- function(plt, breaks, log10) {
  plt <- plt + geom_line() + geom_point() + theme_bw()
  if (log10) {
    plt <- plt + scale_x_log10(breaks = breaks)
  } else {
    plt <- plt + scale_x_continuous(breaks = breaks)
  }
  if (log10) {
    plt <- plt + theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
    plt <- odemodeling:::add_inv_tol_xlab(plt)
  }
  plt + theme(
    panel.grid.minor = element_blank()
  )
}

# Incidence data frame for plotting (SEIR experiment)
get_incidence_quantiles_df <- function(fit) {
  incid <- fit$extract_unflattened(variable = "incidence_gq")
  p <- c(0.1, 0.5, 0.9)
  get_q <- function(x) {
    stats::quantile(x, probs = p)
  }
  aa <- apply(incid, 2, get_q)
  lower <- as.vector(aa[1, ])
  median <- as.vector(aa[2, ])
  upper <- as.vector(aa[3, ])
  t_incid <- t[1:length(t) - 1]
  df <- data.frame(t_incid, lower, median, upper)
  colnames(df) <- c("t", "lower", "median", "upper")
  return(df)
}
