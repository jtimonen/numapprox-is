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

L <- 10 # number of repetitions
times <- matrix(0, L, 2)
ad_calls <- matrix(0, L, 2)
DF <- NULL

for (idx_rep in 1:L) {
  cat(" -------------------- idx_rep =", idx_rep, "--------------------\n")

  # tol = 0.05
  solver1 <- rk4(num_steps = 4)
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

  # tol = 0.02
  solver2 <- rk4(num_steps = 4)
  fit2 <- model$sample(
    t0 = dat$t0_sim, t = dat$t, data = add_data, init = init,
    solver = solver2, step_size = 0.1,
    iter_warmup = ITER,
    iter_sampling = ITER,
    chains = CHAINS
  )
  profs2 <- get_profile_matrix(fit2)
  t2 <- fit2$time()$total
  ad_calls2 <- sum(as.numeric(profs2[6, ]))
  times[idx_rep, ] <- c(t1, t2)
  ad_calls[idx_rep, ] <- c(ad_calls1, ad_calls2)
  pm2 <- rowSums(format_prof_matrix(profs2))

  vals <- c(pm1, pm2)
  names <- as.factor(c(names(pm1), names(pm2)))
  num_steps <- as.factor(rep(c(solver1$num_steps, solver2$num_steps), each = 7))
  df_j <- data.frame(vals, names, num_steps)
  df_j$rep <- as.factor(rep(idx_rep, 7 * 2))
  DF <- rbind(DF, df_j)
}

fn <- file.path(res_dir, "profiling_rk4.rds")
out <- list(times = times, ad_calls = ad_calls, full = DF)
saveRDS(out, fn)
