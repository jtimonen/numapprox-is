# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
library(evmix) # for generalized Pareto distribution density
source("../R/functions.R")
source("tmdd_setup.R")

prior <- tmdd_model(prior_only = TRUE)

# Define simulation parameters
sim_k <- c(0.592, 0.900, 2.212, 0.823, 0.201, 0.024)
sim_sigma <- 0.5
sim_params <- prior$make_params(c(sim_k, sim_sigma))


# Simulate ODE solution
simulate_mp <- function(h) {
  solver <- midpoint(ceiling(10 / h))
  # t_sim <- c(0.1, 0.2, 0.4, 0.6, 1, seq(2, 10, by = 1))
  t_sim <- seq(h, 10, by = h)
  L0_sim <- 10
  t0_sim <- 0
  t_sim <- t_sim
  L0_sim <- L0_sim
  t0_sim <- t0_sim
  sim <- prior$gqs(
    t0 = t0_sim,
    t = t_sim,
    data = list(L0 = L0_sim, D = 3),
    params = sim_params,
    solver = solver
  )
  ynam <- c("y1", "y2", "y3")
  df <- sim$extract_odesol_df(ydim_names = ynam, include_y0 = TRUE)
  df$h <- as.factor(rep(paste0("h = ", h), nrow(df)))
  df
}

# Simulated ODE solution and noisy data
df_sim1 <- simulate_mp(2)
df_sim2 <- simulate_mp(1)
df_sim3 <- simulate_mp(0.5)
df_sim <- rbind(df_sim1, df_sim2, df_sim3)

plt <- ggplot(df_sim, aes(x = t, y = ysol, group = h, color = h)) +
  geom_line() +
  facet_wrap(. ~ ydim) +
  theme_bw() +
  geom_point() +
  ylab("y(t)") +
  xlab("t") +
  theme(legend.title = element_blank())
