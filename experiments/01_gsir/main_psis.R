#!/usr/bin/env Rscript
if (interactive()) {
  idx <- 1
} else {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- args[1]
}
res_dir <- "res_gsir3"

cat("\n------ idx = ", idx, " --------\n", sep = "")
fn_res <- file.path(res_dir, paste0("res_", idx, ".rds"))
fn_run <- file.path(res_dir, paste0("run_", idx, ".rds"))
cat("Results will be loaded from: ", fn_res, "\n", sep = "")
cat("Results will be saved to: ", fn_run, "\n", sep = "")
idx <- as.numeric(idx)

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

# Options
stan_opts <- list(
  sig_figs = 12 # number of significant figures to store in floats
)

# R functions
source("../R/classes.R")
source("../R/functions.R")
source("setup_gsir.R")

# Create experiment setup (need this to recompile Stan models)
solver_args_gen <- list(
  rel_tol = 1e-9,
  abs_tol = 1e-9,
  max_num_steps = 1e9
)
solver <- "rk45"
gpar <- paste("gamma[", c(1:10), "]", sep = "")
param_names <- c("beta", gpar, "phi_inv")
setup <- OdeExperimentSetup$new(
  "gsir", solver, solver_args_gen,
  stan_opts, param_names
)

# Read sampling results
res <- readRDS(file = fn_res)
res$setup$stanmodels <- setup$stanmodels

# Rerun workflow
tol_steps <- 9
run <- rerun_workflow(res, 10, tol_steps, NULL)
# saveRDS(run, file = fn_run)
