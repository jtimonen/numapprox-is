# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
library(ggpubr)


# left plot ---------------------------------------------------------------

# Helper function
time_df <- function(result, ylog) {
  fits <- result$res$fits
  reliab <- result$reliab
  idx <- result$idx
  create_time_comparison_df(fits, reliab, idx, ylog)
}

# Create plot
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
plt_left <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = cols,
    labels = labs
  ) +
  theme(legend.position = c(0.25, 0.7), legend.title = element_blank()) +
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


# right plot --------------------------------------------------------------

out <- results$outputs

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

# Plotting
get_metric_df <- function(output, metric) {
  reliab <- output$reliab
  confs <- output$confs_rel
  if (metric == "max_ratio") {
    plt <- plot_max_ratios_tol(reliab, tols = confs)
  } else {
    plt <- plot_metric_tol(reliab, tols = confs, metric)
  }
  return(plt$data)
}

# Create y-axis label
metric_to_ylabel <- function(metric) {
  if (metric == "max_ratio") {
    str <- paste0("max~~r^{M~','~~M^{'*'}}")
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


# Plot metrics combining different start points
plot_metric_combine <- function(out, metric) {
  cols <- c("#ca0020", "#f4a582", "#92c5de", "#0571b0")
  df <- NULL
  for (j in 1:4) {
    df_j <- get_metric_df(out[[j]], metric)
    df <- rbind(df, df_j)
  }
  aesth <- aes_string(x = "logtol", y = "value", color = "legend")
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_reverse(breaks = unique(round(df$logtol))) +
    xlab("log10(tol)") +
    ylab(metric_to_ylabel(metric)) +
    scale_color_manual(
      values = cols
    ) +
    theme(legend.title = element_blank(), legend.position = c(0.7, 0.7))
  return(plt)
}

plt_A <- plot_metric_combine(out, "mad_odesol")
plt_B <- plot_metric_combine(out, "max_ratio")
plt_C <- plot_metric_combine(out, "pareto_k") +
  geom_hline(yintercept = 0.5, lty = 2)
plt_D <- plot_metric_combine(out, "r_eff")
plt_right <- ggpubr::ggarrange(plt_A, plt_B, plt_C, plt_D)


# combine -----------------------------------------------------------------

plt <- ggpubr::ggarrange(plt_left, plt_right,
  widths = c(0.4, 0.6), labels = "auto"
)

ggsave(plt, file = "figures/tmdd_figure2.pdf", width = 12, height = 4.5)
