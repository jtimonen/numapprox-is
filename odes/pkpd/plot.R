
# Plot data
dat <- readRDS(file='data/dat_sigma_0.1_set_1.rds')
y_hat <- dat$y_hat
y <- dat$y
D <- 5
data_idx <- 1
par(mfrow=c(3,2))
xlim <- c(min(ts), max(ts))

for(d in 1:D){
  ylim <- c(min(y[,d]), max(y[,d]))
  plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[d]))
  lines(ts, y_hat[,d], col=col1, lwd=lwd)
  points(ts, y[,d], col=col2, pch=20)
}



