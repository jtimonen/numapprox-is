library(rstan)
library(gridExtra)
library(tidybayes)
library(readr)
library(dplyr)

# Read data
data_dir <- file.path("disease_transmission_workflow", "data")
fp <- file.path(data_dir, "swiss_agg_data.csv")
df_swiss <- read_csv(fp)

# Cases
cases <- df_swiss$report_dt

# Swiss population
N <- 8.57E6
# initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 <- c(S = s0, I = i0, R = r0)

# times
n_days <- length(cases)
t <- seq(1, n_days, by = 1)
t0 <- 0
t <- t
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)

# Data
data_seir <- list(n_days = n_days, t0 = t0, ts = t, N = N, cases = cases)

# Forcing
date_switch <- "2020-03-13" # date of introduction of control measures
tswitch <- df_swiss %>%
  filter(date < date_switch) %>%
  nrow() + 1 # convert time to number
data_forcing <- list(n_days = n_days, t0 = t0, ts = t, N = N, cases = cases, tswitch = tswitch)
model_forcing_survey <- stan_model("charlesm93_final.stan")

# Survey
date_survey_left <- "2020-05-04"
date_survey_right <- "2020-05-07"
t_survey_start <- df_swiss %>%
  filter(date < date_survey_left) %>%
  nrow() + 1 # convert time to number
t_survey_end <- df_swiss %>%
  filter(date < date_survey_right) %>%
  nrow() + 1 # convert time to number
n_infected_survey <- 83
n_tested_survey <- 775
# add these data to the data given to stan
data_forcing_survey <- c(data_forcing, list(
  t_survey_start = t_survey_start,
  t_survey_end = t_survey_end,
  n_infected_survey = n_infected_survey,
  n_tested_survey = n_tested_survey
))

fit_forcing_survey_max <- sampling(model_forcing_survey,
  data_forcing_survey,
  control = list(max_treedepth = 13, adapt_delta = 0.9, stepsize=0.1),
  iter = 1000,
  seed = 0,
  init = 0,
  refresh = 10
)
