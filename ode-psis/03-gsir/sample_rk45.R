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
  solver = rk45(tol = 1e-5, max_num_steps = 1e4)
)

# Plot
df <- fit$extract_odesol_df(
  ydim_names = c("S", "I", "R", "D"),
  include_y0 = TRUE
)
irow <- which(df$ydim %in% c("I", "D"))
df <- df[irow, ]
draw_inds <- sample.int(size = 100, n = 4000)
irow <- which(df$idx %in% draw_inds)
df <- df[irow, ]
plt <- ggplot(df, aes(x = t, y = ysol, group = idx)) +
  facet_grid(. ~ ydim) +
  geom_line(color = "firebrick", alpha = 0.3)
