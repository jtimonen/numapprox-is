# Setup for the Lotka-Volterra experiments

# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Create model and load data
model <- ode_model_lv()
dat <- load_data_lynxhare()
add_data <- list(
  y_obs_init = dat$y_obs_init,
  y_obs = dat$y_obs,
  D = 2
)

# Define initial params for sampling
init <- list(
  alpha = 1,
  beta = 0.1,
  gamma = 1,
  delta = 0.1,
  y0 = dat$y_obs_init,
  sigma = c(1, 1)
)
init <- rep(list(init), CHAINS) # same for all chains
