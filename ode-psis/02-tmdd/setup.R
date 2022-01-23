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

# Plotting
plot_tmdd <- function(fit, data) {
  y_gen_rvar <- posterior::as_draws_rvars(fit$draws("y_gen"))$y_gen
  df <- create_ribbon_plot_df(y_gen_rvar)
  df$t <- data$t
  cs <- bayesplot::color_scheme_get()
  fill_alpha <- 1.0
  plt_y <- ggplot(df, aes(
    x = t, y = median, ymin = lower1, ymax = upper1
  )) +
    geom_ribbon(alpha = fill_alpha, fill = cs$light_highlight) +
    geom_ribbon(
      alpha = fill_alpha, fill = cs$mid,
      aes(ymin = lower2, ymax = upper2)
    ) +
    geom_line(col = cs$mid_highlight, lwd = 1) +
    ylab("L")
  if (!is.null(data$y)) {
    df <- data.frame(data$t, data$y)
    colnames(df) <- c("t", "y")
    plt_y <- plt_y + geom_point(
      data = df, aes(x = t, y = y),
      inherit.aes = FALSE, pch = 16
    )
  }
  return(plt_y)
}
