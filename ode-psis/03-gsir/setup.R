#!/usr/bin/env Rscript
# Setup GSIR experiment

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Hello
create_data_matrix <- function(x, p) {
  N <- length(x)
  P <- length(p)
  X <- matrix(rep(x, P), nrow = N, byrow = F)
  return(X * p)
}

# Aggregate the groups
aggregate_data <- function(add_data) {
  N <- nrow(add_data$I_data)
  add_data$pop_sizes <- array(sum(add_data$pop_sizes), dim = c(1))
  add_data$G <- 1
  add_data$I_data <- array(rowSums(add_data$I_data), dim = c(N, 1))
  add_data$D_data <- array(rowSums(add_data$D_data), dim = c(N, 1))
  add_data$D <- 4
  add_data$contacts <- matrix(1, 1, 1)
  add_data$I0 <- array(sum(add_data$I0), dim = c(1))
  return(add_data)
}

# Load data
dat <- load_data_lombardia("../../data/lombardia/")
X <- dat$lombardy_data_4may
inc_I <- dat$incidence_cases
inc_D <- dat$incidence_deaths
p_I <- dat$agedistr_cases / sum(dat$agedistr_cases)
p_D <- dat$agedistr_deaths / sum(dat$agedistr_deaths)
pop_sizes <- dat$age_dist * dat$pop_t

# Function
cumulative_column_sums <- function(x){
  y <- array(0, dim=dim(x))
  for(j in 1:ncol(x)){
    y[,j] <- cumsum(x[,j])
  }
  return(y)
}

# Create the actual Stan data
I_mat <- create_data_matrix(inc_I, p_I)
I0 <- I_mat[1, ]
D_mat <- create_data_matrix(inc_D, p_D)
I_data <- round(I_mat[2:nrow(X), ])
D_data <- round(D_mat[2:nrow(X), ])
N <- nrow(I_data)
t0 <- 0
t <- 1:N
G <- length(pop_sizes)
add_data <- list(
  pop_sizes = pop_sizes,
  contacts = dat$contacts,
  G = G,
  I_data = I_data,
  D_data = D_data,
  delta = 0.001,
  I0 = I0,
  D = G * 4,
  deaths_cumulative = cumulative_column_sums(D_data)
)

# Create model
model <- ode_model_gsir()
