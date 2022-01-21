
# Plots
res <- r1
fits <- res$fits
reliab <- res$reliab
plt1_A <- plot_metrics(reliab, tols = res$confs_rel)
plt2_A <- plot_time_comparison_tol(fits, reliab, res$idx)
plt3_A <- plot_time_comparison_tol(fits, reliab, res$idx, TRUE)
diags_A <- get_diags_df(fits) # rhat and reff

# Plots
res <- r2
fits <- res$fits
reliab <- res$reliab
plt1_B <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_B <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_B <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_B <- get_diags_df(fits) # rhat and reff

# Plots
res <- r3
fits <- res$fits
reliab <- res$reliab
plt1_C <- plot_metrics(reliab, num_steps = res$confs_rel)
plt2_C <- plot_time_comparison_ns(fits, reliab, res$idx)
plt3_C <- plot_time_comparison_ns(fits, reliab, res$idx, TRUE)
diags_C <- get_diags_df(fits) # rhat and reff
