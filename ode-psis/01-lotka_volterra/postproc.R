# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)


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

plt <- plot_metrics(reliab, tols = tols_rel)
diags <- get_diags_df(fits)

# Plot times
plot_time_comparison <- function(fits, reliab, idx_ok) {
  tols <- get_tol_vec(fits$solvers)
  tols_rel <- get_tol_vec(reliab$solvers)
  gt <- fits$times$grand_total
  plot(-log10(tols), gt, "o",
    ylab = "time (s)", pch = 16,
    xlab = "T", xaxt = "n"
  )
  grid()
  t_sample <- gt[idx_ok]
  t2 <- reliab$times + t_sample
  lines(-log10(tols_rel), t2, col = "firebrick3")
  points(-log10(tols_rel), t2, col = "firebrick3", pch = 17)
  axis(1, at = -log10(tols), las = 2, labels = tols)
}

plot_time_comparison(fits, reliab, idx)
