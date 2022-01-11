# Requirements
library(cmdstanr)
library(posterior)
library(bayesplot)
library(checkmate)
library(loo)
library(stats)
library(outbreaks)
library(scales)
library(ggplot2)
library(ggdist)
library(R6)
library(ggpubr)
library(tidyr)
source("../R/classes.R")
source("../R/functions.R")
source("setup_gsir.R")

R <- 60
ntol <- 11
sim_times <- NULL # tuning time
sampling_times <- NULL # sampling time
full_times <- NULL # init + sampling time
k_hats <- NULL
maes <- NULL
for (idx in 1:60) {
  cat("idx=", idx, "\n", sep = "")
  tryCatch(
    {
      fn_res <- file.path("res_gsir3", paste0("res_", idx, ".rds"))
      fn_run <- file.path("res_gsir3", paste0("run_", idx, ".rds"))
      res <- readRDS(file = fn_res)
      run <- readRDS(file = fn_run)
      tim <- res$run$post_fit$time()
      tim_sum_chains <- sum(tim$chains$total)
      times_sum_chains <- rowSums(res$tps$total)
      sim_times <- rbind(sim_times, run$tuning$metrics$time)
      sampling_times <- rbind(sampling_times, c(tim_sum_chains, times_sum_chains))
      full_times <- rbind(full_times, c(tim$total, res$tps$grand_total))
      sim_tols <- 1 / run$tuning$metrics$inv_tol
      k_hats <- rbind(k_hats, run$tuning$metrics$k_hat)
      maes <- rbind(maes, run$tuning$metrics$mae)
    },
    error = function(e) {
      cat(" - Could not open, idx = \n", idx, ".\n", sep = "")
    }
  )
}
sampling_tols <- 10^c(-3, -4, -6, -8, -10, -12)
colnames(sim_times) <- sim_tols
colnames(k_hats) <- sim_tols
colnames(maes) <- sim_tols
colnames(sampling_times) <- sampling_tols
colnames(full_times) <- sampling_tols
diffs <- full_times - sampling_times
mns <- rep(c(1e3, 1e4, 1e5), each = 20)
pltA <- plot_timing(sim_tols, sim_times, mns) + ggtitle("Log prob computation time (no gradients)")
pltB <- plot_timing(sampling_tols, full_times, mns, unit = "hours") + ggtitle("Sampling time")
pltC <- plot_timing(sampling_tols, diffs, mns) + ggtitle("sampler.init_stepsize() time")
plt1 <- ggarrange(pltA, pltB, pltC, nrow = 3, labels = c("A", "B", "C"))

# Speedup plot
su <- compute_speedup(sim_tols, sampling_tols, sim_times, full_times)
su <- su[, 2:ncol(su)]
plt2 <- plot_timing(NULL, su, mns = mns) +
  ylab("(t_sample_low + t_prob_eval_high) / t_sample_high") + ggtitle("Relative time, shown number is median")

# Reliability
pltA <- plot_metric(k_hats, mns) + ylab("Pareto-k")
pltB <- plot_metric(maes, mns) + ylab("Max absolute error")
plt3 <- ggarrange(pltA, pltB, nrow = 2, labels = c("A", "B"))

# Save plots
ggsave(filename = "figs/timing.pdf", plot = plt1, width = 9.3, height = 10)
ggsave(filename = "figs/relative_time.pdf", plot = plt2, width = 12, height = 5.5)
ggsave(filename = "figs/reliability.pdf", plot = plt3, width = 8.5, height = 10)
