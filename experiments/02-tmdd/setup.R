# Setup for the TMDD experiments

# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

# Create model and load data
dat <- readRDS("simulated_data.rds")
model <- ode_model_tmdd()
add_data <- list(L0 = dat$L0, P_obs = dat$P_obs, D = 3)
init <- 0
