#!/usr/bin/env Rscript

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Setup data and model
source("setup.R")

# Sampling
fit <- model$sample(
  t0 = t0, t = t, data = add_data, init = init,
  solver = rk45(tol = 1e-4, max_num_steps = 1e4),
  step_size = step_size, iter_warmup = 300, iter_sampling = 200, chains = 1
)

max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
print(max_rhat)

# Plot
df <- fit$extract_odesol_df(
  ydim_names = c("S", "E", "I", "R"),
  include_y0 = TRUE
)
irow <- which(df$ydim %in% c("I", "D"))
df <- df[irow, ]
draw_inds <- sample.int(size = 100, n = 4000)
irow <- which(df$idx %in% draw_inds)
df <- df[irow, ]
plt <- ggplot(df, aes(x = t, y = ysol, group = idx)) +
  facet_grid(. ~ ydim, scales = "free_y") +
  geom_line(color = "firebrick", alpha = 0.3)

dead <- add_data$deaths_cumulative
df_dead <- data.frame(t = t, y = dead, ydim = rep("D", length(dead)))
df_dead$ydim <- as.factor(df_dead$ydim)
plt <- plt + geom_point(
  data = df_dead,
  aes(x = t, y = y),
  inherit.aes = FALSE
)
