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
  max_rhat <- max(fit$summary()[, "rhat"])
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
