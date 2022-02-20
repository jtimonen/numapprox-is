# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)
library(scales)

# Load results
res_dir <- c("results_bdf")
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)

# Plotting
plot_results <- function(res, ylog = TRUE) {
  fits <- res$res$fits
  reliab <- res$reliab
  list(
    metrics = plot_metrics(reliab, tols = res$confs_rel),
    times = plot_time_comparison_tol(fits, reliab, res$idx, ylog),
    diags = get_diags_df(fits) # rhat and reff
  )
}

# Create plots
p1 <- plot_results(results$outputs[[1]])
p2 <- plot_results(results$outputs[[2]])
p3 <- plot_results(results$outputs[[3]])
p4 <- plot_results(results$outputs[[4]])
