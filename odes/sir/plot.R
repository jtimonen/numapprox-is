# Plot options
col1     <- 'gray70'
col2     <- 'gray10'
lwd      <- 2

# Plot data
dat <- readRDS(file='data/dat_set_1.rds')
y_hat <- dat$y_hat
y <- dat$y
D <- 20
data_idx <- 1
par(mfrow=c(1,1))
par(mar = c(2,2,2,2) + 0.1)
ts <- dat$ts
xlim <- c(min(ts), max(ts))

for(d in 2:2){
  ylim <- c(min(y), max(y))
  plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[d]))
  lines(ts, y_hat[,d], col=col1, lwd=lwd)
  points(ts, y, col=col2, pch=20)
}
