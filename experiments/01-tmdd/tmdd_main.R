#!/usr/bin/env Rscript
# TMDD experiment

# Run all TMDD experiments
# - requires that tmdd_data.rds exists

# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(loo)
library(evmix) # for generalized Pareto distribution density
source("../R/functions.R")

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results"
odemodeling:::create_dir_if_not_exist(res_dir)
source("tmdd_setup.R") # defines data and model

# MCMC sampling, takes a couple of hours
source("tmdd_mcmc.R")

# PSIS
source("tmdd_reliability.R")

# Create result figures
source("figures/create_figure_simulation.R")
source("figures/create_figure_mainresults.R")
source("figures/create_figure_psis.R")

