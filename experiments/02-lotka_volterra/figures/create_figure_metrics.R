# Load results
dirs <- c("rk45", "rk4", "midpoint")
par_dirs <- file.path(res_dir, dirs)
results <- list()
for (j in 1:3) {
  par_dir <- par_dirs[j]
  fp <- file.path(par_dir, "reliability.rds")
  results[[j]] <- readRDS(file = fp)
}
names(results) <- dirs


# RK45 metrics --------------------------------------------------------------

# Plot metrics
plot_metric_rk45 <- function(out, metric) {
  df <- get_metric_df_tol(out, metric)
  aesth <- aes_string(x = "logtol", y = "value", color = "legend")
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    xlab("log10(tol)") +
    ylab(metric_to_ylabel(metric)) +
    theme(legend.title = element_blank(), legend.position = c(0.55, 0.4))
  return(plt)
}

plt_A <- plot_metric_rk45(out, "mad_odesol")
plt_B <- plot_metric_rk45(out, "max_ratio")
plt_C <- plot_metric_rk45(out, "pareto_k")
plt_D <- plot_metric_rk45(out, "r_eff")
plt_rk45_metrics <- ggpubr::ggarrange(plt_A, plt_B, plt_C, plt_D)



# RK4 and midpoint metrics -------------------------------------------------

# Plot metrics for midpoint and RK4
plot_metric_numsteps <- function(out, metric, color) {
  df <- get_metric_df_ns(out, metric)

  aesth <- aes_string(x = "num_steps", y = "value", color = "legend")
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab("K") +
    scale_x_continuous(breaks = seq(2, 30, by = 2)) +
    ylab(metric_to_ylabel(metric)) +
    scale_color_manual(values = color) +
    theme(legend.title = element_blank(), legend.position = c(0.55, 0.4))
  return(plt)
}

# Combine
plot_4_metrics <- function(out, color) {
  plt_A <- plot_metric_numsteps(out, "mad_odesol", color)
  plt_B <- plot_metric_numsteps(out, "max_ratio", color)
  plt_C <- plot_metric_numsteps(out, "pareto_k", color)
  plt_D <- plot_metric_numsteps(out, "r_eff", color)
  ggpubr::ggarrange(plt_A, plt_B, plt_C, plt_D, nrow = 1)
}

p_rk4 <- plot_4_metrics(results$rk4, "black")
p_mp <- plot_4_metrics(results$midpoint, "black")

# Combine -----------------------------------------------------------------
