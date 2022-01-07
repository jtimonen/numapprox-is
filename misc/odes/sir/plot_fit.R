# Plot options
col1 <- "gray70"
col2 <- "gray10"
col3 <- "firebrick"
col4 <- "orange"
col5 <- "steelblue"
lwd <- 2

dat <- readRDS(file = "data/dat_sigma_1_set_1.rds")
fit <- readRDS(file = "data/yhat.rds")

y_hat <- dat$y_hat
y <- dat$y
D <- 20
data_idx <- 1
par(mfrow = c(1, 2))
par(mar = c(2, 2, 2, 2) + 0.1)
ts <- dat$ts
xlim <- c(min(ts), max(ts))

ylim <- c(min(y), max(y))

i_smp <- 23
for (d in 1:2) {

  # Plot data
  plot(0, 0, ylim = ylim, xlim = xlim, pch = NA, xlab = "t", ylab = expression(y[, d]))
  lines(ts, y_hat[, d], col = col1, lwd = lwd)
  points(ts, y[, d], col = col2, pch = 20)

  # Plot fit
  for (j in i_smp) {
    yi <- fit$inf[j, , d]
    yr <- fit$ref[j, , d]
    lines(ts, yr, col = col3, lwd = 2)
    lines(ts, yi, col = col4, lwd = 2)
    points(ts, yi, col = col4, pch = 4)
  }
  points(ts, y[, d], col = col2, pch = 20)
}



# Another plot
# plot(0, 0, xlim=xlim, ylim = c(0, 1000))
# for(j in 1:100){
#  lines(ts, fit$ref[j,,1], col=col4)
#  lines(ts, fit$ref[j,,2], col=col3)
#  lines(ts, fit$ref[j,,3], col=col5)
# }
