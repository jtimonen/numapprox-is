#!/usr/bin/env Rscript
# Setup SIRD experiment

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Compile model
model <- ode_model_seir()

# Load data
add_data <- load_data_switzerland("../../data/switzerland/")
