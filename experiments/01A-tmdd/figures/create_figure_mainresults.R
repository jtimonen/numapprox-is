# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)

# Create timing plot
ylog <- FALSE
df <- NULL
lab0 <- expression(time[MCMC]^{
  BDF(tol)
})
labs <- list()
labs[[1]] <- lab0
for (j in 1:4) {
  out <- results$outputs[[j]]
  tol_bdf <- out$confs[out$idx]
  df_j <- time_df(out, ylog)
  df_j$logtol <- log10(1 / df_j$inv_tol)
  if (j > 1) {
    df_j <- df_j[which(df_j$procedure != "high"), ]
  }
  df_j$procedure <- as.character(df_j$procedure)
  df_j$procedure[which(df_j$procedure != "high")] <- paste0("low", j)
  df <- rbind(df, df_j)
  str <- paste0("time[MCMC]^{BDF(", tol_bdf, ")} + time[PSIS]^{BDF(tol)}")
  labs[[j + 1]] <- parse(text = str)
}
df$procedure <- as.factor(df$procedure)

# Plot
n_yticks <- 8
cols <- c("#010101", "#ca0020", "#f4a582", "#92c5de", "#0571b0")
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_timing <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = cols,
    labels = labs
  ) +
  theme(legend.position = c(0.2, 0.6), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df$logtol))) +
  xlab("log10(tol)") +
  ylab("time (s)") +
  scale_y_continuous(
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10^x,
      n = n_yticks
    ),
    labels = trans_format("log10", math_format(10^.x))
  )

# Metrics plots
out <- results$outputs

# Plot metrics
plot_metric_bdf <- function(out, metric, color = "#4daf4a") {
  df <- get_metric_df_tol(out, metric)
  aesth <- aes_string(x = "logtol", y = "value", color = "legend")
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    scale_color_manual(values = color) +
    xlab("log10(tol)") +
    ylab(metric_to_ylabel(metric)) +
    theme(legend.title = element_blank(), legend.position = c(0.55, 0.4),
          panel.grid.minor = element_blank())
  return(plt)
}

plot_bdf <- function(out, color = "#4daf4a") {
  plt_A <- plot_metric_bdf(out, "mad_odesol", color) +
    theme(legend.position = "none")
  plt_B <- plot_metric_bdf(out, "max_log_ratio", color) +
    theme(legend.position = "none")
  plt_C <- plot_metric_bdf(out, "pareto_k", color) +
    theme(legend.position = "none")
  plt_D <- plot_metric_bdf(out, "r_eff", color)
  ggpubr::ggarrange(plt_A, plt_B, plt_C, plt_D, nrow = 1)
}

p1 <- plot_bdf(out[[1]], color = cols[2])
p2 <- plot_bdf(out[[2]], color = cols[3])
p3 <- plot_bdf(out[[3]], color = cols[4])
p4 <- plot_bdf(out[[4]], color = cols[5])
plt_metrics <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 1)

ggsave(plt_timing, file = "figures/tmdd_timing.pdf", width = 6, height = 4.5)
ggsave(plt_metrics, file = "figures/tmdd_metrics.pdf", width = 12, height = 7.5)
