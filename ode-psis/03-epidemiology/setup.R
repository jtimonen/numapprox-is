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
dat <- load_data_switzerland("../../data/switzerland/")

# Create actual Stan data
add_data <- dat[c(
  "n_days", "cases", "tswitch", "t_survey_start",
  "t_survey_end", "n_infected_survey", "n_tested_survey"
)]

add_data$pop_size <- dat$N
add_data$D <- 4
add_data$n_days_m1 <- add_data$n_days - 1
t <- dat$ts
t0 <- dat$t0

# Other
init <- 0
step_size <- 0.1
