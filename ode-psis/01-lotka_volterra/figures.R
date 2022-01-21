# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Load results
dirs <- c("results_rk45", "results_rk4", "results_midpoint")
results <- list()
for (j in 1:3) {
  res_dir <- dirs[j]
  fp <- file.path(res_dir, "reliability.rds")
  results[[j]] <- readRDS(file = fp)
}
names(results) <- dirs

# Plotting
plot_results <- function(res, ylog = TRUE) {
  fits <- res$res$fits
  reliab <- res$reliab
  is_adaptive <- is(fits$solvers[[1]], "AdaptiveOdeSolver")
  reliab <- res$reliab
  if (is_adaptive) {
    plt_A <- plot_metrics(reliab, tols = res$confs_rel)
    plt_B <- plot_time_comparison_tol(fits, reliab, res$idx, ylog)
  } else {
    plt_A <- plot_metrics(reliab, num_steps = res$confs_rel)
    plt_B <- plot_time_comparison_ns(fits, reliab, res$idx, ylog)
  }
  list(
    metrics = plt_A,
    times = plt_B,
    diags = get_diags_df(fits) # rhat and reff
  )
}

# Create plots
p1 <- plot_results(results[[1]])
p2 <- plot_results(results[[2]])
p3 <- plot_results(results[[3]])
