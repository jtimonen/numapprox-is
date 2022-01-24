# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Load results
dirs <- c("results_bdf")
results <- list()
inds <- c(1)
for (j in 1:1) {
  index <- inds[j]
  res_dir <- dirs[j]
  fn <- paste0("reliability_idx", index, ".rds")
  fp <- file.path(res_dir, fn)
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

# BDF --------------------------------------------------------------------

tol_bdf <- results[[1]]$confs[results[[1]]$idx]

df1 <- time_df(results[[1]], ylog)
df1$logtol <- log10(1 / df1$inv_tol)
lab1 <- expression(time[MCMC]^{
  BDF(tol)
})
lab2 <- expression(time[MCMC]^
  {
    BDF(0.06)
  } + time[IS]^{
    BDF(tol)
  })

labs <- c(lab1, lab2)
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_A <- ggplot(df1, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_discrete(labels = labs) +
  theme(legend.position = c(0.7, 0.25), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df1$logtol))) +
  xlab("log10(tol)")
if (ylog) {
  plt_A <- plt_A + ylab("log(time)")
}

# Combine
plt <- plt_A
ggsave(plt, filename = "figures/times.pdf", width = 4, height = 3.5)
