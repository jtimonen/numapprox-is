# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
source("../R/functions.R")
source("tmdd_setup.R")

# cmdstanr::set_cmdstan_path("/Users/juhotimonen/Work/research/stan/cmdstan")
ITER <- 2000
CHAINS <- 4
model <- tmdd_model(profile_oderhs = TRUE)
res_dir <- "results_profiling" # results won't be saved
odemodeling:::create_dir_if_not_exist(res_dir)

#  Helper function
get_profile_matrix <- function(fit) {
  p <- sapply(fit$cmdstanr_fit$profiles(), function(x) x)
  as.matrix(p[3:nrow(p), ])
}

#  Helper function
format_prof_matrix <- function(x) {
  p <- matrix(as.numeric(x), 7)
  rownames(p) <- rownames(x)
  p
}

# tol = 0.05
solver1 <- bdf(tol = 1e-6, max_num_steps = 1e5)
fit1 <- model$sample(
  t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
  solver = solver1, step_size = 0.1,
  iter_warmup = ITER,
  iter_sampling = ITER,
  chains = CHAINS
)
profs1 <- get_profile_matrix(fit1)
t1 <- fit1$time()$total
ad_calls1 <- sum(as.numeric(profs1[6, ]))
pm1 <- rowSums(format_prof_matrix(profs1))
time_per_grad <- t1 / ad_calls1

# Run gq
gq1 <- fit1$gqs()
t1_gq <- gq1$time()$total
time_per_nograd <- t1_gq / gq1$ndraws()

fn <- file.path(res_dir, "profiling_grad_no_grad.rds")
out <- list(times = times, ad_calls = ad_calls, full = DF)
# saveRDS(out, fn)

cat("grad:", time_per_grad, "\n")
cat("no_grad:", time_per_nograd, "\n")
