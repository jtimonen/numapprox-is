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
add_data <- aggregate_data(add_data)

# Sampling
fit <- model$sample(
  t0 = t0, t = t, data = add_data, init = 0,
  solver = rk45(tol = 1e-4, max_num_steps = 1e5)
)
