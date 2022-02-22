#!/usr/bin/env Rscript
# Lotka-Volterra experiment

# Requirements
library(odemodeling)
library(posterior)
library(scales)
library(ggplot2)
source("../R/functions.R")

# Setup
args <- commandArgs(trailingOnly = TRUE)
ITER <- 2000
CHAINS <- 4
res_dir <- "results"
odemodeling:::create_dir_if_not_exist(res_dir)
source("lv_setup.R") # defines data and model

# Sampling
source("lv_mcmc_midpoint.R")
source("lv_mcmc_rk4.R")
source("lv_mcmc_rk45.R")

# PSIS
source("lv_reliability.R")

# Figures
source("figures/create_figure1.R")
