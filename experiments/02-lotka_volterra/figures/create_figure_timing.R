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


# RK45 ----------------------------------------------------------------------

# Create plot for RK45
out <- results$rk45
lab0 <- expression(time[MCMC]^{
  RK45(tol)
})

tol_rk45 <- out$confs[out$idx]
df <- time_df(out, FALSE)
df$logtol <- log10(1 / df$inv_tol)
df$procedure <- as.character(df$procedure)
df$procedure[which(df$procedure != "high")] <- "low"
str <- paste0("time[MCMC]^{RK45(", tol_rk45, ")} + time[PSIS]^{RK45(tol)}")
labs <- list(lab0, parse(text = str))
df$procedure <- as.factor(df$procedure)

# Plot
ymax <- 10^3.6
ymin <- 10^1.6
n_yticks <- 8
cols_rk45 <- c("black", "#4daf4a")
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_rk45_time <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = ymax, alpha = 0) +
  geom_hline(yintercept = ymin, alpha = 0) +
  scale_color_manual(
    values = cols_rk45,
    labels = labs
  ) +
  theme(legend.position = c(0.25, 0.7), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df$logtol))) +
  xlab("log10(tol)") +
  ylab("time (s)") +
  scale_y_continuous(
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10^x, n = n_yticks),
    labels = trans_format("log10", math_format(10^.x))
  )


# RK4 and midpoint --------------------------------------------------------
cols_ns <- c("orange", "#984ea3", "#e41a1c", "#377eb8")
ns_rk4 <- results[[2]]$confs[results[[2]]$idx]
ns_mp <- results[[3]]$confs[results[[3]]$idx]
df2 <- time_df(results[[2]], FALSE)
df3 <- time_df(results[[3]], FALSE)
df <- rbind(df2, df3)
method <- c(rep("rk4", nrow(df2)), rep("midpoint", nrow(df3)))
df$method <- as.factor(method)
df$group <- paste0(df$procedure, " (", df$method, ")")

lab1a <- expression(time[MCMC]^{
  RK4(K)
})
lab1b <- expression(time[MCMC]^{
  midpoint(K)
})
lab2a <- expression(time[MCMC]^
  {
    RK4(2)
  } + time[PSIS]^{
    RK4(K)
  })
lab2b <- expression(time[MCMC]^
  {
    midpoint(3)
  } + time[PSIS]^{
    midpoint(K)
  })

labs <- c(lab1b, lab1a, lab2b, lab2a)
aesth <- aes(x = num_steps, y = time, group = group, color = group)
plt_ns_time <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = ymax, alpha = 0) +
  geom_hline(yintercept = ymin, alpha = 0) +
  theme(legend.position = c(0.7, 0.45), legend.title = element_blank()) +
  scale_color_manual(
    values = cols_ns,
    labels = labs,
  ) +
  xlab("Number of steps K") +
  ylab("time (s)") +
  scale_x_continuous(breaks = seq(2, 30, by = 2)) +
  scale_y_continuous(
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10^x, n = n_yticks),
    labels = trans_format("log10", math_format(10^.x))
  )


# combine -----------------------------------------------------------------

plt <- ggarrange(plt_rk45_time, plt_ns_time, labels = "auto")
ggsave(plt, filename = "figures/lv_timing.pdf", width = 7.75, height = 3.43)
