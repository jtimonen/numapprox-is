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

# Plot
ynam <- c("x1", "x2", "x3")
plt_sim <- sim$plot_odesol(
  ydim_names = ynam, include_y0 = TRUE
)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep(ynam[3], length(t_sim)))
plt_sim <- plt_sim +
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab("Concentration") + xlab("t")


# inference results -------------------------------------------------------

# Load results
res_dir <- c("results_bdf")
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
out <- results$outputs[[3]]
fits <- out$res$fits
fit <- readRDS(fits$files[1])
fit$plot_odesol_dist(include_y0 = TRUE)
