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

# Create better plots
odemodeling:::create_dir_if_not_exist("figures")

# Helper function
time_df <- function(result, ylog) {
  fits <- result$res$fits
  reliab <- result$reliab
  idx <- result$idx
  create_time_comparison_df(fits, reliab, idx, ylog)
}
ylog <- TRUE

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
  xlab("log10(tol)") +
  ylim(3.5, 8)
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
  scale_x_continuous(breaks = seq(2, 30, by = 2)) +
  ylim(3.5, 8)
if (ylog) {
  plt_B <- plt_B + ylab("log(time)")
}

# Combine
plt <- ggarrange(plt_A, plt_B, labels = "auto")
ggsave(plt, filename = "figures/times.pdf", width = 7.75, height = 3.43)
