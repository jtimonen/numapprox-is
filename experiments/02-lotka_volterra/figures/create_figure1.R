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
library(ggpubr)


# left plot ---------------------------------------------------------------

# Helper function
time_df <- function(result) {
  fits <- result$res$fits
  reliab <- result$reliab
  idx <- result$idx
  create_time_comparison_df(fits, reliab, idx, FALSE)
}

# Create plot for RK45
out <- results$rk45
lab0 <- expression(time[MCMC]^{
  RK45(tol)
})

tol_rk45 <- out$confs[out$idx]
df <- time_df(out)
df$logtol <- log10(1 / df$inv_tol)
df$procedure <- as.character(df$procedure)
df$procedure[which(df$procedure != "high")] <- "low"
str <- paste0("time[MCMC]^{RK45(", tol_rk45, ")} + time[PSIS]^{RK45(tol)}")
labs <- list(lab0, parse(text = str))
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


# RK45 --------------------------------------------------------------------

tol_rk45 <- results[[1]]$confs[results[[1]]$idx]
df1 <- time_df(results[[1]], ylog)
df1$logtol <- log10(1 / df1$inv_tol)
lab1 <- expression(time[MCMC]^{
  RK45(tol)
})
lab2 <- expression(time[MCMC]^
  {
    RK45(0.001)
  } + time[IS]^{
    RK45(tol)
  })

labs <- c(lab1, lab2)
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_A <- ggplot(df1, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_discrete(labels = labs) +
  theme(legend.position = c(0.7, 0.45), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df1$logtol))) +
  xlab("log10(tol)")
if (ylog) {
  plt_A <- plt_A + ylab("log(time)")
}


# RK4 and midpoint --------------------------------------------------------

ns_rk4 <- results[[2]]$confs[results[[2]]$idx]
ns_mp <- results[[3]]$confs[results[[3]]$idx]
df2 <- time_df(results[[2]], ylog)
df3 <- time_df(results[[3]], ylog)
df <- rbind(df2, df3)
method <- c(rep("rk4", nrow(df2)), rep("midpoint", nrow(df3)))
df$method <- as.factor(method)
df$group <- paste0(df$procedure, " (", df$method, ")")

lab1a <- expression(time[MCMC]^{
  RK4(M)
})
lab1b <- expression(time[MCMC]^{
  midpoint(M)
})
lab2a <- expression(time[MCMC]^
  {
    RK4(2)
  } + time[IS]^{
    RK4(M)
  })
lab2b <- expression(time[MCMC]^
  {
    midpoint(3)
  } + time[IS]^{
    midpoint(M)
  })

labs <- c(lab1b, lab1a, lab2b, lab2a)
aesth <- aes(x = num_steps, y = time, group = group, color = group)
plt_B <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_brewer(type = "div", palette = 5, labels = labs) +
  xlab("Number of steps M") +
  theme(legend.position = c(0.7, 0.45), legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(2, 30, by = 2))
if (ylog) {
  plt_B <- plt_B + ylab("log(time)")
}

# Combine
plt <- ggarrange(plt_A, plt_B, labels = "auto")
ggsave(plt, filename = "figures/times.pdf", width = 7.75, height = 3.43)
