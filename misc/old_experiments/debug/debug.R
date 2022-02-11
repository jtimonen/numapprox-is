library(cmdstanr)
set_cmdstan_path("C:\\Users/Juho/Work/Research/STAN/cmdstan/")

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
  chains = 2,
  save_warmup = TRUE,
  step_size = 0.1 # 0.01
)
diag <- fit$sampler_diagnostics(inc_warmup = TRUE)
diag[, , c("treedepth__", "n_leapfrog__", "stepsize__")]

print(fit$time())
print(fit$profiles())
