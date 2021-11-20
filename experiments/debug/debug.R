
library(cmdstanr)

dat <- readRDS(file = "gsir_data.rds")
model <- cmdstan_model("gsir_post.stan")
fit <- model$sample(
  data = dat,
  sig_figs = 12,
  seed = 30785100,
  refresh = 1,
  show_messages = TRUE,
  init = 0,
  iter_warmup = 1,
  iter_sampling = 1,
  chains = 1
)
