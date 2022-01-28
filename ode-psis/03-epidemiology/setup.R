#!/usr/bin/env Rscript
# Setup SIRD experiment

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Load data
dat <- load_data_lombardia("../../data/lombardia/")
X <- dat$lombardy_data_4may
I0 <- dat$incidence_cases[1]
N <- length(dat$incidence_cases)
dead <- cumsum(dat$incidence_deaths)[2:N]
N <- length(dead)

# Create the actual Stan data
t0 <- 0
t <- 1:N
add_data <- list(
  pop_size = dat$pop_t,
  contacts = dat$contacts,
  delta = 0.001,
  I0 = I0,
  D = 4,
  deaths_cumulative = dead
)

# Create model
model <- ode_model_sird()
