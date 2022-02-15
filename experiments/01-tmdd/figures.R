# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)


# simulation --------------------------------------------------------------

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
P_dat <- simdat$P_obs

# Simulation ODE solution and noisy data
ynam <- c("x1", "x2", "x3")
df_sim <- sim$extract_odesol_df(ydim_names = ynam, include_y0 = TRUE)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep(ynam[3], length(t_sim)))

# Load results
res_dir <- c("results_bdf")
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
out <- results$outputs[[3]]
fits <- out$res$fits

# Load two fits and plot their ODE solution distribution
fit_low <- readRDS(fits$files[3]) # mcmc fit with tol=0.04
fit_high <- readRDS(fits$files[16]) # mcmc fit with tol=1e-12
probs <- c(0.05, 0.5, 0.95)
qlow <- fit_low$extract_odesol_df_dist(
  p = probs, ydim_names = ynam, include_y0 = TRUE
)
qhigh <- fit_high$extract_odesol_df_dist(
  p = probs, ydim_names = ynam, include_y0 = TRUE
)
df_dist <- rbind(qlow, qhigh)
tol_low <- fit_low$solver$abs_tol
tol_high <- fit_high$solver$abs_tol
df_dist$tol <- as.factor(rep(c(tol_low, tol_high), each = nrow(qlow)))
colnames(df_dist) <- c("t", "ydim", "lower", "median", "upper", "tol")
plt <- ggplot(df_dist, aes(
  x = t, y = median, group = tol, color = tol,
  fill = tol, ymin = lower, ymax = upper, lty = tol
)) +
  geom_line() +
  facet_wrap(. ~ ydim) +
  geom_ribbon(alpha = 0.1) +
  theme_bw() +
  geom_line(data = df_sim, aes(x = t, y = ysol), inherit.aes = FALSE) +
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab("Concentration") +
  xlab("t")

# ggsave(plt, file = "figures/tmdd_simdata.pdf", width=7.9, height=6.3)
