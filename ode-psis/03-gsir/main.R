#!/usr/bin/env Rscript
# TMDD experiment
args <- commandArgs(trailingOnly = TRUE)

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")

library(odemodeling)
library(posterior)

# Setup
idx <- setup_experiment_index(args)
fp <- setup_experiment_paths(idx)

dat <- load_data_lombardia("../../data/lombardy/")
