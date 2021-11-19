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
source("classes.R")
source("functions.R")
source("setup_gsir.R")

R <- 60
ntol <- 11
ttimes <- NULL
ptimes <- NULL
for (idx in 1:60) {
  cat("idx=", idx, "\n", sep = "")
  tryCatch(
    {
      fn <- file.path("res_gsir", paste0("res_", idx, ".rds"))
      res <- readRDS(file = fn)
      ttimes <- rbind(ttimes, res$run$tuning$metrics$time)
      ptimes <- rbind(ptimes, rowSums(res$tps$total))
      tols <- 1 / res$run$tuning$metrics$inv_tol
    },
    error = function(e) {
      cat(" - Could not open!\n")
    }
  )
}
colnames(ttimes) <- tols
colnames(ptimes) <- tols
