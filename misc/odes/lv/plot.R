#!/usr/bin/env Rscript
library(pracma) # for errorbar plot

# Plot options
col1 <- "gray70"
col2 <- "gray10"
lwd <- 2


# Create figure
fig_name <- paste0("data.pdf")
pdf(fig_name, width = 12, height = 14)
par(mfrow = c(5, 4))

for (data_idx in 1:10) {
  # Load data
  sigma <- 0.5
  data_idx <- 3
  fn <- paste0("data/dat_set_", data_idx, ".rds")
  dat <- readRDS(file = fn)
  y <- dat$y
  y_hat <- dat$y_hat
  ts <- dat$ts

  # Plot data
  ylim <- c(min(y), max(y))
  xlim <- c(min(ts), max(ts))
  plot(0, 0, ylim = ylim, xlim = xlim, pch = NA, xlab = "t", ylab = expression(y[1]))
  lines(ts, y_hat[, 1], col = col1, lwd = lwd)
  points(ts, y[, 1], col = col2, pch = 20)

  plot(0, 0, ylim = ylim, xlim = xlim, pch = NA, xlab = "t", ylab = expression(y[2]))
  lines(ts, y_hat[, 2], col = col1, lwd = lwd)
  points(ts, y[, 2], col = col2, pch = 20)
}


# Close graphics device
dev.off()
