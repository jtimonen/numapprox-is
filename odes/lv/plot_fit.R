# Plot options
col1 <- "gray70"
col2 <- "gray10"
col3 <- "firebrick"
col4 <- "orange"
col5 <- "steelblue"
lwd <- 2

# Plot data
dat <- readRDS(file = "data/dat_set_1.rds")
y_hat <- dat$y_hat
y <- dat$y
D <- 20
data_idx <- 1
par(mfrow = c(1, 1))
par(mar = c(2, 2, 2, 2) + 0.1)
ts <- dat$ts
xlim <- c(min(ts), max(ts))

ylim <- c(min(y), max(y))
plot(0, 0, ylim = ylim, xlim = xlim, pch = NA, xlab = "t", ylab = expression(y[2]))
lines(ts, y_hat[, 2], col = col1, lwd = lwd)
points(ts, y, col = col2, pch = 20)

fit <- readRDS(file = "data/yhat.rds")

# Plot fit
for (j in 1:200) {
  yf <- fit$ref[j, , 2]
  lines(ts, yf, col = col3)
}
points(ts, y, col = col2, pch = 20)

# Another plot
# plot(0, 0, xlim=xlim, ylim = c(0, 1000))
# for(j in 1:100){
#  lines(ts, fit$ref[j,,1], col=col4)
#  lines(ts, fit$ref[j,,2], col=col3)
#  lines(ts, fit$ref[j,,3], col=col5)
# }
