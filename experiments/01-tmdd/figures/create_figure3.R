# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)

# Plotting
plot_results <- function(res, ylog = TRUE) {
  fits <- res$res$fits
  reliab <- res$reliab
  list(
    metric = plot_metric_tol(reliab, tols = res$confs_rel, "pareto_k"),
    max_ratio = plot_max_ratios_tol(reliab, tols = res$confs_rel),
    diags = get_diags_df(fits) # rhat and reff
  )
}

# Create plots
p1 <- plot_results(results$outputs[[1]])
p2 <- plot_results(results$outputs[[2]])
p3 <- plot_results(results$outputs[[3]])
p4 <- plot_results(results$outputs[[4]])
