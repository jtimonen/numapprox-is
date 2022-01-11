#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Setup -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ITER <- 20
CHAINS <- 2

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Setup
idx <- setup_experiment_index(args)
fp <- setup_experiment_paths(idx)

# Create model and load data
model <- ode_model_lv()
dat <- load_data_lynxhare()
add_data <- list(
  y_obs_init = dat$y_obs_init,
  y_obs = dat$y_obs,
  D = 2
)

# Define initial params for sampling
init <- list(
  alpha = 1,
  beta = 0.1,
  gamma = 1,
  delta = 0.1,
  y0 = dat$y_obs_init,
  sigma = c(1, 1)
)
init <- rep(list(init), CHAINS) # same for all chains

# SAMPLING ----------------------------------------------------------
# solver_fit <- rk45(abs_tol = 1e-5, rel_tol = 1e-5, max_num_steps = 1e5)
tols <- c(0.01, 0.005, 10^(-3:-12))
solvers <- rk45_list(tols = tols, max_num_steps = 1e4)
fits <- model$sample_manyconf(
  t0 = dat$t0, t = dat$t, data = add_data, init = init,
  solvers = solvers, step_size = 0.1, savedir = fp$res_dir,
  iter_warmup = ITER, iter_sampling = ITER, chains = CHAINS
)

# Load first fist
fit <- readRDS(file = fits$files[1])

gq <- fit$gqs(solver = rk45(tol = 1e-12))
is <- psis(fit, gq)
print(is$diagnostics)

# Run workflow
idx <- 2
sampled <- load_fit(post, idx)
L <- length(tols)
tols_val <- tols[(idx + 1):L]
run <- validate_fit(setup, sampled, tols_val, max_num_steps)

# Save result
tune <- run$tuning
tp <- plot_tuning(tune)
ggsave("tmdd_tuning.pdf", width = 9.67, height = 6.17)

# Plot times
t <- post$grand_total[idx:L]
tols <- post$tols[idx:L]
plot(-log10(tols), t, "o",
  ylab = "time (s)", pch = 16, ylim = c(0, 320),
  xlab = "T", xaxt = "n"
)
grid()
t_sample <- t[1]
t2 <- tune$time + t_sample
tols2 <- 1 / tune$inv_tol
lines(-log10(tols2), t2, col = "firebrick3")
points(-log10(tols2), t2, col = "firebrick3", pch = 17)
legend(2, 200, c("HMC-NUTS using tol=T", "HMC-NUTS using tol=0.05 + PSIS with tol=T"),
  lty = c(1, 1), col = c("black", "firebrick3"),
  pch = c(16, 17)
)
axis(1, at = -log10(tols), las = 2, labels = tols)
