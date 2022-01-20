# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)
library(ggplot2)


# RK45 --------------------------------------------------------------------

# Read in sampling results
res_rk45 <- readRDS("results_rk45/sampling.rds")
fits <- res_rk45$fits
tols <- get_tol_vec(fits$solvers)

# Load fit
idx <- 3 # 1 and 2 fail
inds_rel <- (1 + idx):length(tols)
fit <- load_fit(file = fits$files[idx])

# Run reliability check
tols_rel <- tols[inds_rel]
rel_solvers <- rk45_list(tols = tols_rel, max_num_steps = 1e9)
reliab <- fit$reliability(
  solvers = rel_solvers, force = TRUE,
  savedir = "results_rk45"
)

# Plot reliability metrics
plt <- plot_metrics(reliab, tols = tols_rel)
diags <- get_diags_df(fits) # rhat and reff

# Plot times
plt2 <- plot_time_comparison_tol(fits, reliab, idx)
plt3 <- plot_time_comparison_tol(fits, reliab, idx, TRUE)


# RK4 ---------------------------------------------------------------------

# Read in sampling results
res_rk4 <- readRDS("results_rk4/sampling.rds")
fits <- res_rk4$fits
ns <- get_num_steps_vec(fits$solvers)

# Load fit
idx <- 1 # 1 and 2 fail
inds_rel <- (1 + idx):length(ns)
fit <- load_fit(file = fits$files[idx])

# Run reliability check
ns_rel <- ns[inds_rel]
rel_solvers <- rk4_list(num_steps = ns_rel)
reliab <- fit$reliability(
  solvers = rel_solvers, force = TRUE,
  savedir = "results_rk4"
)

# Plot reliability metrics
plt <- plot_metrics(reliab, num_steps = ns_rel)
diags <- get_diags_df(fits) # rhat and reff

# Plot times
plt2 <- plot_time_comparison_ns(fits, reliab, idx)
plt3 <- plot_time_comparison_ns(fits, reliab, idx, TRUE)
