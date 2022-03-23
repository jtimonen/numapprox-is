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

L <- 10 # number of repetitions
times <- matrix(0, L, 2)
ad_calls <- matrix(0, L, 2)
for (idx_rep in 1:L) {

  # tol = 0.05
  solver1 <- bdf(tol = 0.05, max_num_steps = 1e5)
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

  # tol = 0.02
  solver2 <- bdf(tol = 0.02, max_num_steps = 1e5)
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
}

fn <- file.path(res_dir, "profiling.rds")
out <- list(times = times, ad_calls = ad_calls)
saveRDS(out, fn)

tol1 <- solver1$abs_tol
tol2 <- solver2$abs_tol
tols <- as.factor(rep(c(tol1, tol2), each = L))
df <- data.frame(as.vector(times), as.vector(ad_calls), tols)
colnames(df) <- c("time", "ad_calls", "tol")
plt <- ggplot(df, aes(x = ad_calls, y = time, group = tol, color = tol)) +
  geom_smooth(
    method = "lm", aes(x = ad_calls, y = time), inherit.aes = FALSE,
    color = "black", se = FALSE, lwd = 0.5
  ) +
  theme_bw() +
  xlab("Number of ODE RHS calls with autodiff") +
  ylab("time (s)") +
  geom_point() +
  scale_color_brewer(type = "qual", palette = 6)

ggsave(plt, filename = "figures/profiling.pdf", width = 5.7, height = 4)
