# Load simulated data
simdat <- readRDS(file = "tmdd_data.rds")

# Create model and simulation solver
prior <- tmdd_model(prior_only = TRUE)

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
P_dat <- simdat$P_obs

# Simulation ODE solution and noisy data
ynam <- c("y1", "y2", "y3")
df_sim <- sim$extract_odesol_df(ydim_names = ynam, include_y0 = TRUE)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep(ynam[3], length(t_sim)))

# Plot simulated data and ODE solution
plt_sim <- ggplot(data = df_sim, aes(x = t, y = ysol), inherit.aes = FALSE) +
  geom_line() +
  facet_wrap(. ~ ydim) +
  theme_bw() +
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab(paste("epsilon =", sim$solver$abs_tol))

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
idx_out <- 4
out <- results$outputs[[idx_out]]
idx_low <- results$start_inds[idx_out]
fits <- out$res$fits
rel_files <- out$reliab$files
gq_low <- out$reliab$base
gq_high <- readRDS(file = rel_files[length(rel_files)])

# Absolute errors
maes <- abs(gq_low$extract_odesol() - gq_high$extract_odesol())
maes <- apply(maes, 1, max)
a <- sort(maes, index.return = T, decreasing = TRUE)
inds <- a$ix[1:10]

# Plot low-accuracy solutions
# inds <- sample.int(gq_low$ndraws(), size = 20, replace = FALSE)
plt_low <- gq_low$plot_odesol(draw_inds = inds) +
  theme_bw() + ylab(paste("epsilon =", gq_low$solver$abs_tol))
plt_high <- gq_high$plot_odesol(draw_inds = inds) +
  theme_bw() + ylab(paste("epsilon =", gq_high$solver$abs_tol))

# Combine
plt <- ggarrange(plt_sim, plt_low, plt_high, nrow = 3, labels = "auto")
ggsave(plt, file = "figures/tmdd_simulation.pdf", width = 8.1, height = 7.1)
