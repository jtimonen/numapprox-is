library(ggplot2)

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

# Plot max likelihood ratios
plot_max_ratios_tol <- function(rel, tols) {
  tol_base <- rel$base$solver$abs_tol
  leg <- paste0("tol_star=", tol_base)
  L <- length(tols)
  legend <- rep(leg, L)

  # Compute all likelihood ratios
  log_ratios <- all_log_ratios(rel)
  max_ratios <- apply(exp(log_ratios), 2, max)

  # Create data frame
  df <- data.frame(log10(tols), max_ratios, legend)
  solver_name <- toupper(rel$base$solver$name)
  labs <- paste0(
    "M = ", solver_name, "(", tol_base,
    "),  M* = ", solver_name, "(tol)"
  )
  df$legend <- factor(legend, labels = labs)
  colnames(df) <- c("logtol", "value", "legend")

  # Create y label
  str <- paste0("max~~r^{M~','~~M^{'*'}}")
  ylabel <- parse(text = str)

  # Plot
  aesth <- aes(x = logtol, y = value, color = legend)
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    xlab("log10(tol)") +
    ylab(ylabel)
  return(plt)
}

# Plot any other metric
plot_metric_tol <- function(rel, tols, metric) {
  tol_base <- rel$base$solver$abs_tol
  leg <- paste0("tol_star=", tol_base)
  L <- length(tols)
  legend <- rep(leg, L)

  # Compute all likelihood ratios
  values <- rel$metrics[, metric]

  # Create data frame
  df <- data.frame(log10(tols), values, legend)
  solver_name <- toupper(rel$base$solver$name)
  labs <- paste0(
    "M = ", solver_name, "(", tol_base,
    "),  M* = ", solver_name, "(tol)"
  )
  df$legend <- factor(legend, labels = labs)
  colnames(df) <- c("logtol", "value", "legend")

  aesth <- aes(x = logtol, y = value, color = legend)
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    xlab("log10(tol)") +
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
