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
df2 <- time_df(results[[2]], ylog)
df3 <- time_df(results[[3]], ylog)
df <- rbind(df2, df3)
method <- c(rep("rk4", nrow(df2)), rep("midpoint", nrow(df3)))
df$method <- as.factor(method)
df$group <- paste0(df$procedure, " (", df$method, ")")

plt <- ggplot(df, aes(x = num_steps, y = time, group = group, 
                      color = group)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(type = "div", palette = 5, 
                     labels = c("a", "b", "c", "d")) +
  theme_bw()
if (ylog) {
  plt <- plt + ylab("log(time)")
}

lab1 <- expression(tau[1])
lab2 <- expression(tau[2])
labs <- c(lab1, lab2)