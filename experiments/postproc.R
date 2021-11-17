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

idx <- 1
fn <- file.path("res_gsir", paste0("res_", idx, ".rds"))
res <- readRDS(file = fn)

tols <- 1/res$run$tuning$metrics$inv_tol
plot_timing(tols, res$tps$total)
