# Setup for the TMDD experiments

# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)

dat <- readRDS("simulated_data.rds")

model <- ode_model_tmdd()

data <- list(L0 = 10, D = 3, t = dat$t, P_obs = dat$P_obs)

# Data
setup_standata_tmdd <- function() {
  t <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0)
  N <- length(t)
  L0 <- 14.8
  data_list <- list(
    N = N,
    t = t,
    L0 = L0,
    D = 3,
    y = rep(1.0, N) # dummy
  )
  return(data_list)
}

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
