# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Load simulated data
simdat <- readRDS(file = "simulated_data.rds")

# Create model and simulation solver
prior <- ode_model_tmdd(prior_only = TRUE)

# Define simulation parameters
sim_k <- simdat$sim_k
sim_sigma <- simdat$sim_sigma
sim_params <- prior$make_params(c(sim_k, sim_sigma))

# Simulate ODE solution
t_sim <- simdat$t
L0_sim <- simdat$L0_sim
t0_sim <- simdat$t0_sim
sim <- prior$gqs(
  t0 = t0_sim,
  t = t_sim,
  data = list(L0 = L0_sim, D = 3),
  params = sim_params,
  solver = simdat$solver_sim
)
y_sol <- sim$extract_odesol()
P_dat <- simdat$P_obs

# Plot
ynam <- c("x1", "x2", "x3")
plt_sim <- sim$plot_odesol(
  ydim_names = ynam, include_y0 = TRUE
)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep(ynam[3], length(t_sim)))
plt_sim <- plt_sim +
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab("Concentration") + xlab("t")

# Load results
res_dir <- c("results_bdf")
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)


# Plotting
plot_results <- function(res, ylog = TRUE) {
  fits <- res$res$fits
  reliab <- res$reliab
  list(
    metrics = plot_metrics(reliab, tols = res$confs_rel),
    times = plot_time_comparison_tol(fits, reliab, res$idx, ylog),
    diags = get_diags_df(fits) # rhat and reff
  )
}

# Create plots
p1 <- plot_results(results$outputs[[1]])
p2 <- plot_results(results$outputs[[2]])
p3 <- plot_results(results$outputs[[3]])
p4 <- plot_results(results$outputs[[4]])

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
  str <- paste0("time[MCMC]^{BDF(", tol_bdf, ")} + time[IS]^{BDF(tol)}")
  labs[[j + 1]] <- parse(text = str)
}
df$procedure <- as.factor(df$procedure)

# Plot
cols <- c("#010101", "#ca0020", "#f4a582", "#92c5de", "#0571b0")
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_A <- ggplot(df, aesth) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = cols,
    labels = labs
  ) +
  theme(legend.position = c(0.2, 0.75), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df$logtol))) +
  xlab("log10(tol)")
if (ylog) {
  plt_A <- plt_A + ylab("log(time)")
}

# Combine
plt <- plt_A
# ggsave(plt, filename = "figures/times.pdf", width = 4, height = 3.5)
