
library(cmdstanr)

dat <- readRDS(file = "gsir_data.rds")
model <- cmdstan_model("gsir_post.stan")
#dat$abs_tol <- 1e-6
#dat$rel_tol <- 1e-6
fit <- model$sample(
  data = dat,
  sig_figs = 12,
  seed = 30785100,
  refresh = 1,
  show_messages = TRUE,
  init = 0,
  iter_warmup = 1,
  iter_sampling = 1,
  chains = 1,
  save_warmup = TRUE,
  save_latent_dynamics=TRUE
)
diag <- fit$sampler_diagnostics(inc_warmup = TRUE)
diag[,,c("treedepth__", "n_leapfrog__", "stepsize__")]
