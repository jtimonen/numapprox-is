library(ggplot2)
library(ggpubr)

# Helper function
create_log10_epsilon_label <- function() {
  parse(text = paste0("log10(", expression(epsilon), ")"))
}

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

# Plotting helper
get_metric_df_tol <- function(output, metric) {
  reliab <- output$reliab
  confs <- output$confs_rel
  if (metric == "max_log_ratio") {
    plt <- plot_max_ratios_tol(reliab, tols = confs)
  } else {
    plt <- plot_metric_tol(reliab, tols = confs, metric)
  }
  return(plt$data)
}

# Plotting helper
get_metric_df_ns <- function(output, metric) {
  reliab <- output$reliab
  confs <- output$confs_rel
  if (metric == "max_log_ratio") {
    plt <- plot_max_ratios_ns(reliab, num_steps = confs)
  } else {
    plt <- plot_metric_ns(reliab, num_steps = confs, metric)
  }
  return(plt$data)
}

# Create y-axis label
metric_to_ylabel <- function(metric) {
  if (metric == "max_log_ratio") {
    str <- paste0("max~~log~r^{M~','~~M^{'*'}}")
    ylabel <- parse(text = str)
  } else if (metric == "mad_odesol") {
    str <- paste0("max~~'|'~y^{M}-y^{M^{'*'}}~'|'")
    ylabel <- parse(text = str)
  } else if (metric == "pareto_k") {
    ylabel <- "Pareto-k"
  } else if (metric == "r_eff") {
    ylabel <- "Relative efficiency"
  } else {
    ylabel <- metric
  }
  return(ylabel)
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

# Max likelihood ratios
all_log_ratios <- function(rel) {
  J <- length(rel$files)
  base <- rel$base
  a <- list()
  for (j in 1:J) {
    gq <- readRDS(rel$files[j])
    a[[j]] <- as.vector(log_ratios(base, gq))
  }
  no_op <- function(x) x
  sapply(a, no_op)
}

# Plot max likelihood ratios (tol on x axis)
plot_max_ratios_tol <- function(rel, tols) {
  tol_base <- rel$base$solver$abs_tol
  solver_name <- toupper(rel$base$solver$name)

  # Compute all likelihood ratios
  log_ratios <- all_log_ratios(rel)
  max_log_ratios <- apply(log_ratios, 2, max)
  values <- max_log_ratios

  # Create y label
  str <- paste0("max~~log~r^{M~','~~M^{'*'}}")
  metric <- parse(text = str)

  plot_metric_tol_impl(tol_base, tols, values, metric, solver_name)
}

# Plot any other metric (tol on x axis)
plot_metric_tol <- function(rel, tols, metric) {
  tol_base <- rel$base$solver$abs_tol
  solver_name <- toupper(rel$base$solver$name)
  values <- rel$metrics[, metric]
  plot_metric_tol_impl(tol_base, tols, values, metric, solver_name)
}

create_labels_tol <- function(solver_name, tol_base) {
  paste0(
    "M = ", solver_name, "(", tol_base,
    "),\nM* = ", solver_name, "(epsilon)"
  )
}

plot_metric_tol_impl <- function(tol_base, tols, values, metric, solver_name) {
  leg <- paste0("tol_base=", tol_base)
  L <- length(tols)
  legend <- rep(leg, L)

  # Create data frame
  df <- data.frame(log10(tols), values, legend)
  labs <- create_labels_tol(solver_name, tol_base)
  df$legend <- factor(legend, labels = labs)
  colnames(df) <- c("logtol", "value", "legend")

  aesth <- aes(x = logtol, y = value, color = legend)
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    xlab(create_log10_epsilon_label()) +
    ylab(metric)
  return(plt)
}


# Plot max likelihood ratios (num_steps on x axis)
plot_max_ratios_ns <- function(rel, num_steps) {

  # Compute all likelihood ratios
  log_ratios <- all_log_ratios(rel)
  max_log_ratios <- apply(log_ratios, 2, max)
  values <- max_log_ratios

  # Create y label
  str <- paste0("max~~log~r^{M~','~~M^{'*'}}")
  metric <- parse(text = str)

  ns_base <- rel$base$solver$num_steps
  solver_name <- toupper(rel$base$solver$name)
  plot_metric_ns_impl(ns_base, num_steps, values, metric, solver_name)
}

# Plot any other metric (num steps on x axis)
plot_metric_ns <- function(rel, num_steps, metric) {
  ns_base <- rel$base$solver$num_steps
  values <- rel$metrics[, metric]
  solver_name <- toupper(rel$base$solver$name)
  plot_metric_ns_impl(ns_base, num_steps, values, metric, solver_name)
}

create_labels_ns <- function(solver_name, ns_base) {
  paste0("M = ", solver_name, "(", ns_base, "),\nM* = ", solver_name, "(K)")
}

plot_metric_ns_impl <- function(ns_base, num_steps, values, metric, solver_name) {
  leg <- paste0("ns_base=", ns_base)
  L <- length(num_steps)
  legend <- rep(leg, L)

  # Create data frame
  df <- data.frame(num_steps, values, legend)
  labs <- create_labels_ns(solver_name, ns_base)
  df$legend <- factor(legend, labels = labs)
  colnames(df) <- c("num_steps", "value", "legend")

  aesth <- aes(x = num_steps, y = value, color = legend)
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("K") +
    ylab(metric)
  return(plt)
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

# Helper function
time_df <- function(result, ylog) {
  fits <- result$res$fits
  reliab <- result$reliab
  idx <- result$idx
  create_time_comparison_df(fits, reliab, idx, ylog)
}
