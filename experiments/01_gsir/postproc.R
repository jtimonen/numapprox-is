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
    },
    error = function(e) {
      cat(" - Could not open, idx = \n", idx, ".\n", sep = "")
    }
  )
}
sampling_tols <- 10^c(-3, -4, -6, -8, -10, -12)
colnames(sim_times) <- sim_tols
colnames(sampling_times) <- sampling_tols
colnames(full_times) <- sampling_tols
diffs <- full_times - sampling_times
mns <- rep(c(1e3, 1e4, 1e5), each = 20)
plt1 <- plot_timing(sim_tols, sim_times, mns) + ggtitle("Tuning time")
plt2 <- plot_timing(sampling_tols, full_times, mns, hours = TRUE) + ggtitle("Sampling time")
plt3 <- plot_timing(sampling_tols, diffs, mns) + ggtitle("sampler.init_stepsize() time")
plt <- ggarrange(plt1, plt2, plt3, nrow = 3, labels = c("A", "B", "C"))

# Speedup plot
asdd

# reliability plot
asd
