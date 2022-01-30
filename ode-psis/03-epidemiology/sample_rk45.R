#!/usr/bin/env Rscript

# R functions and requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)

# Setup data and model
source("setup.R")
ITER <- 1000
print(model$stanmodel)

# Sampling
fit <- model$sample(
  t0 = t0, t = t, data = add_data, init = 0,
  solver = rk45(rel_tol = 1e-4, abs_tol = 1e-4, max_num_steps = 1e2),
  step_size = step_size, iter_warmup = ITER, iter_sampling = ITER, chains = 4,
  adapt_delta = 0.9, max_treedepth = 13
)

max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
print(max_rhat)
print(fit$summary)


# Plot incidence
df_dat <- data.frame(dat$cases, dat$ts)
colnames(df_dat) <- c("cases", "t")
df <- get_incidence_quantiles_df(fit)
map <- aes(x = t, y = cases)
plt <- ggplot(df, aes(x = t, y = median, ymin = lower, ymax = upper)) +
  geom_line(color = "firebrick") +
  geom_ribbon(fill = "firebrick2", alpha = 0.6) +
  geom_point(data = df_dat, mapping = map, inherit.aes = FALSE) +
  ylab("Number of reported cases")
