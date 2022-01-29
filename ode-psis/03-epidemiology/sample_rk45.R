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
  t0 = t0, t = t, data = add_data, init = 1,
  solver = midpoint(), # rk45(tol = 1e-2, max_num_steps = 1e4),
  step_size = step_size, iter_warmup = 300, iter_sampling = 200, chains = 4
)

max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
print(max_rhat)

# Plot incidence
df_dat <- data.frame(dat$cases, dat$ts)
colnames(df_dat) <- c("cases", "t")
df <- get_incidence_quantiles_df(fit)
plt <- ggplot(df, aes(x = t, y = median, ymin = lower, ymax = upper)) +
  geom_line(color = "firebrick") +
  geom_ribbon(fill = "firebrick2", alpha = 0.6) +
  geom_point(data = df_dat, mapping = aes(x = t, y = cases), inherit.aes = FALSE) +
  ylab("Incidence")
