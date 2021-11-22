library(cmdstanr)
set_cmdstan_path("~/Work/Research/cmdstan")

dat <- readRDS(file = "gsir_data.rds")
model <- cmdstan_model("gsir_post.stan")
# dat$abs_tol <- 1e-6
# dat$rel_tol <- 1e-6
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
  save_latent_dynamics = TRUE,
  step_size = 1 # 0.01
)
diag <- fit$sampler_diagnostics(inc_warmup = TRUE)
diag[, , c("treedepth__", "n_leapfrog__", "stepsize__")]

print(fit$time())
print(fit$profiles())
