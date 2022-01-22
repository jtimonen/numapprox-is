#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results_bdf"
odemodeling:::create_dir_if_not_exist(res_dir)
source("setup.R")


# SAMPLING ----------------------------------------------------------
max_num_steps <- 1e6
tols <- c(
  0.1, 0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
  1e-11, 1e-12
)
post <- setup$sample_posterior_many(res_dir, idx, tols, max_num_steps, chains = 4)

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


# Run sampling
tols <- c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 10^(-5:-12))
solvers <- rk45_list(tols = tols, max_num_steps = 1e5)
fits <- model$sample_manyconf(
  t0 = dat$t0, t = dat$t, data = add_data, init = init,
  solvers = solvers, step_size = 0.1, savedir = res_dir,
  iter_warmup = ITER, iter_sampling = ITER, chains = CHAINS
)

# Save results
results <- list(
  fits = fits,
  session_info = sessionInfo()
)
fp <- file.path(res_dir, "sampling.rds")
saveRDS(results, file = fp)
