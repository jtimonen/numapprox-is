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

# Plots
res <- results[[1]]
fits <- res$res$fits
reliab <- res$reliab
plt1_A <- plot_metrics(reliab, tols = res$confs_rel)
plt2_A <- plot_time_comparison_tol(fits, reliab, res$idx)
plt3_A <- plot_time_comparison_tol(fits, reliab, res$idx, TRUE)
diags_A <- get_diags_df(fits) # rhat and reff

# Plots
res <- results[[2]]
fits <- res$res$fits
reliab <- res$reliab
plt1_B <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_B <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_B <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_B <- get_diags_df(fits) # rhat and reff

# Plots
res <- results[[3]]
fits <- res$res$fits
reliab <- res$reliab
plt1_C <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_C <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_C <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_C <- get_diags_df(fits) # rhat and reff
