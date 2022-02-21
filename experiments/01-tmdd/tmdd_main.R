#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Run all TMDD experiments
# - requires that simulated_data.rds exists

# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(ggplot2)

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results"
odemodeling:::create_dir_if_not_exist(res_dir)
source("tmdd_setup.R") # defines data and model

# MCMC sampling, takes long
source("tmdd_mcmc.R")

# PSIS, is quick
source("tmdd_reliability.R")

# Create result figures
source("figures/create_figure1.R")
source("figures/create_figure2.R")
# source("figures/create_figure3.R")
