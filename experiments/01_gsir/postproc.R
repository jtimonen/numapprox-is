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
ttimes <- NULL # tuning time
stimes <- NULL # sampling time
gtimes <- NULL # init + sampling time
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
      ttimes <- rbind(ttimes, run$tuning$metrics$time)
      stimes <- rbind(stimes, c(tim_sum_chains, times_sum_chains))
      gtimes <- rbind(gtimes, c(tim$total, res$tps$grand_total))
      ttols <- 1 / run$tuning$metrics$inv_tol
    },
    error = function(e) {
      cat(" - Could not open!\n")
    }
  )
}
stols <- 10^c(-3, -4, -6, -8, -10, -12)
colnames(ttimes) <- ttols
colnames(stimes) <- stols
colnames(gtimes) <- stols
diffs <- gtimes - stimes
mns <- rep(c(1e3, 1e4, 1e5), each = 20)
plt1 <- plot_timing(ttols, ttimes, mns)
plt2 <- plot_timing(stols, stimes, mns)
plt3 <- plot_timing(stols, gtimes, mns)
plt4 <- plot_timing(stols, diffs, mns)

ggarrange(plotlist = plt1, nrow = 1)
# Plot
diff <- dtimes[, 1] - dtimes[, 2]
plot(diff,
  pch = 16, col = "firebrick",
  ylab = "total - sum(chains$total)"
)
grid()
