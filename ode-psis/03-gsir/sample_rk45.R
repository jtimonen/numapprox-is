#!/usr/bin/env Rscript
# Setup GSIR experiment

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Setup data and model
source("setup.R")

# Sampling
fit <- model$sample(t0=t0, t=t, data = add_data)