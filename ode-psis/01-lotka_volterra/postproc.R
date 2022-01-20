# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)
library(ggplot2)

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
