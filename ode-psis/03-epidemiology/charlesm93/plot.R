library(rstan)
library(gridExtra)
library(tidybayes)
library(readr)
library(dplyr)

# Load
out <- readRDS("fit_forcing_survey_max.rds")
fit_forcing_survey_max <- out$fit

# Data
n_days <- out$dat$n_days
cases <- out$dat$cases
c_posterior <- "firebrick"

# Plot posterior
smr_pred <- cbind(as.data.frame(summary(fit_forcing_survey_max, pars = "pred_cases", probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))$summary), t = 1:(n_days - 1), cases = cases[1:length(cases) - 1])
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
plt_post <- ggplot(smr_pred, mapping = aes(x = t)) +
  # geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_dark, ) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = c_posterior, alpha = 0.35) +
  # geom_ribbon(aes(ymin = X10., ymax = X90.), fill = c_light) +
  geom_line(mapping = aes(x = t, y = X50.), color = c_posterior) +
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Incidence")


# Plot Reff
tswitch <- out$dat$tswitch
plt_R0 <- fit_forcing_survey_max %>%
  spread_draws(Reff[n_days]) %>%
  group_by(n_days) %>%
  summarise(R0_mean = mean(Reff), R09 = quantile(Reff, 0.95), R01 = quantile(Reff, 0.05)) %>%
  ggplot() +
  geom_ribbon(aes(x = n_days, ymin = R01, ymax = R09), fill = c_posterior, alpha = 0.35) +
  geom_line(mapping = aes(n_days, R0_mean), color = c_posterior) +
  geom_vline(aes(xintercept = tswitch))
