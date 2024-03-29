# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
library(evmix) # for generalized Pareto distribution density
source("../R/functions.R")
source("tmdd_setup.R")

# Create model and simulation solver
solver_sim <- bdf(
  rel_tol = 1e-15,
  abs_tol = 1e-15,
  max_num_steps = 1e9
)
prior <- tmdd_model(prior_only = TRUE)

# Define simulation parameters
sim_k <- c(0.592, 0.900, 2.212, 0.823, 0.201, 0.024)
sim_sigma <- 0.5
sim_params <- prior$make_params(c(sim_k, sim_sigma))

# Simulate ODE solution
t_sim <- c(0.1, 0.2, 0.4, 0.6, 1, seq(2, 10, by = 1))
L0_sim <- 10
t0_sim <- 0

# Alternative
# sim_sigma <- 0.1
# sim_params <- prior$make_params(c(sim_k, sim_sigma))
# t_sim <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, seq(2, 10, by = 1))
# L0_sim <- 3
# t0_sim <- 0

sim <- prior$gqs(
  t0 = t0_sim,
  t = t_sim,
  data = list(L0 = L0_sim, D = 3),
  params = sim_params,
  solver = solver_sim
)
y_sol <- sim$extract_odesol()
P_sol <- y_sol[1, , 3]

# Add noise and save data
SEED <- 324
set.seed(SEED)
P_dat <- P_sol + rnorm(P_sol, sd = sim_sigma)
dat <- list(
  t = t_sim, P_sol = P_sol, P_obs = P_dat, sim_k = sim_k,
  sim_sigma = sim_sigma, solver_sim = solver_sim, L0_sim = L0_sim,
  t0_sim = t0_sim, seed = SEED
)
saveRDS(dat, file = "tmdd_data.rds")


# plotting ----------------------------------------------------------------

# Load simulated data
simdat <- readRDS(file = "tmdd_data.rds")

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


# Simulated ODE solution and noisy data
ynam <- c("y1", "y2", "y3")
df_sim <- sim$extract_odesol_df(ydim_names = ynam, include_y0 = TRUE)
df_dat <- data.frame(t = t_sim, y = P_dat, ydim = rep(ynam[3], length(t_sim)))

plt <- ggplot(df_sim, aes(x = t, y = ysol)) +
  geom_line() +
  facet_wrap(. ~ ydim) +
  theme_bw() +
  geom_point(data = df_dat, aes(x = t, y = y), inherit.aes = FALSE) +
  ylab("Concentration") +
  xlab("t") +
  theme(legend.title = element_blank())

ggsave(plt, file = "figures/tmdd_simdata.pdf", width = 8.2, height = 3)
