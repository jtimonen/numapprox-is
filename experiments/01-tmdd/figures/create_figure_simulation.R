# Load simulated data
simdat <- readRDS(file = "tmdd_data.rds")

# Create model and simulation solver
prior <- tmdd_model(prior_only = TRUE)

# Define simulation parameters
sim_k <- simdat$sim_k
sim_sigma <- simdat$sim_sigma
sim_params <- prior$make_params(c(sim_k, sim_sigma))

# Simulate ODE solution
t_data <- simdat$t
t_dense <- seq(0.01, max(t_data), by = 0.01)
L0_sim <- simdat$L0_sim
t0_sim <- simdat$t0_sim
sim <- prior$gqs(
  t0 = t0_sim,
  t = t_data,
  data = list(L0 = L0_sim, D = 3),
  params = sim_params,
  solver = simdat$solver_sim
)
sim_dense <- prior$gqs(
  t0 = t0_sim,
  t = t_dense,
  data = list(L0 = L0_sim, D = 3),
  params = sim_params,
  solver = simdat$solver_sim
)
P_dat <- simdat$P_obs

# Plotting simulation ODE solution and noisy data
plot_simul_with_data <- function(sim) {
  ynam <- c("y1", "y2", "y3")
  df_sim <- sim$extract_odesol_df(ydim_names = ynam, include_y0 = TRUE)
  df_dat <- data.frame(t = t_data, y = P_dat, ydim = rep(ynam[3], length(t_data)))

  # Plot simulated data and ODE solution
  ggplot(data = df_sim, aes(x = t, y = ysol), inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(. ~ ydim) +
    theme_bw() +
    geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
    ylab(paste("epsilon =", sim$solver$abs_tol))
}
plt_sim <- plot_simul_with_data(sim)
plt_sim_dense <- plot_simul_with_data(sim_dense)

# Get ODE solution at given time point
get_timepoint <- function(sim, t, draw_inds = NULL) {
  ynam <- c("y1", "y2", "y3")
  df_sim <- sim$extract_odesol_df(
    ydim_names = ynam, include_y0 = TRUE,
    draw_inds = draw_inds
  )
  df_sim[df_sim$t == t, ]
}

# Get absolute difference of simulations at given time point
simul_diff <- function(sim1, sim2, t, draw_inds = NULL) {
  df1 <- get_timepoint(sim1, t, draw_inds = draw_inds)
  df2 <- get_timepoint(sim2, t, draw_inds = draw_inds)
  abs(df1$ysol - df2$ysol)
}

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
idx_out <- 4
out <- results$outputs[[idx_out]]
idx_low <- results$start_inds[idx_out]
fits <- out$res$fits
rel_files <- out$reliab$files

# Load results using low accuracy
gq_low <- out$reliab$base
fit_low <- readRDS(fits$files[out$idx])
fit_low$model$reinit()

# Data input for simulation
data_sim <- list(
  P_obs = rep(1, length(t_data)),
  L0 = L0_sim,
  D = 3
)

# Data input for dense simulation
data_sim_dense <- list(
  P_obs = rep(1, length(t_dense)),
  L0 = L0_sim,
  D = 3
)

# Confirm that resimulation gives same result
sim_low <- fit_low$gqs(
  t0 = t0_sim,
  t = t_data,
  data = data_sim
)
mae1 <- max_abs_odesol_diff(sim_low, gq_low)
message("mae_low_vs_resim = ", mae1, "\n")

# Dense low-accuracy simulation
sim_low_dense <- fit_low$gqs(
  t0 = t0_sim,
  t = t_dense,
  data = data_sim_dense
)

# Load results using high accuracy
gq_high <- readRDS(file = rel_files[length(rel_files)])

# Confirm that resimulation gives same result
sim_high <- fit_low$gqs(
  t0 = t0_sim,
  t = t_data,
  data = data_sim,
  solver = gq_high$solver
)
mae2 <- max_abs_odesol_diff(sim_high, gq_high)
message("mae_high_vs_resim = ", mae2, "\n")

# Dense high-accuracy simulation
sim_high_dense <- fit_low$gqs(
  t0 = t0_sim,
  t = t_dense,
  solver = gq_high$solver,
  data = data_sim_dense
)

# Absolute errors between low and high
maes <- abs(gq_low$extract_odesol() - gq_high$extract_odesol())
maes <- apply(maes, 1, max)
a <- sort(maes, index.return = T, decreasing = TRUE)
inds <- a$ix[1:10] # determine draw inds

# Plot solutions
plot_solutions <- function(gq, inds) {
  gq$plot_odesol(draw_inds = inds, include_y0 = TRUE) +
    theme_bw() + ylab(paste("epsilon =", gq$solver$abs_tol))
}
# inds <- sample.int(gq_low$ndraws(), size = 20, replace = FALSE)
plt_low <- plot_solutions(gq_low, inds)
plt_high <- plot_solutions(gq_high, inds)
plt_low_dense <- plot_solutions(sim_low_dense, inds)
plt_high_dense <- plot_solutions(sim_high_dense, inds)

# Plot at data points
plt <- ggarrange(plt_sim, plt_low, plt_high, nrow = 3, labels = "auto")
ggsave(plt, file = "figures/tmdd_simulation.pdf", width = 8.1, height = 7.1)

# Plot at dense
plt_dense <- ggarrange(plt_sim_dense, plt_low_dense, plt_high_dense, nrow = 3, labels = "auto")
ggsave(plt_dense, file = "figures/tmdd_simulation_dense.pdf", width = 8.1, height = 7.1)


# Differences of dense and non-dense simulation at t = 10
diff_low <- simul_diff(sim_low, sim_low_dense, t = 10, draw_inds = inds)
diff_high <- simul_diff(sim_high, sim_high_dense, t = 10, draw_inds = inds)
cat("max_abs_diff_high = ", max(diff_high), "\n", file = "figures/diff.txt")
cat("max_abs_diff_low = ", max(diff_low), "\n",
  file = "figures/diff.txt",
  append = TRUE
)
