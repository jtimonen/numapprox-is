
# Plot options
col1     <- 'gray70'
col2     <- 'gray10'
lwd      <- 2

# Load data
sigma    <- 0.1
data_idx <- 1
fn       <- paste0('data/dat_sigma_', sigma, '_set_', data_idx, '.rds')
dat      <- readRDS(file=fn)
y        <- dat$y
y_hat    <- dat$y_hat
ts       <- dat$ts

# Plot data
jpeg('dat.jpg')
ylim <- c(min(y), max(y))
xlim <- c(min(ts), max(ts))
par(mfrow=c(2,1))
plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[1]))
lines(ts, y_hat[,1], col=col1, lwd=lwd)
points(ts, y[,1], col=col2, pch=20)

plot(0, 0, ylim=ylim, xlim=xlim, pch=NA, xlab='t', ylab=expression(y[2]))
lines(ts, y_hat[,2], col=col1, lwd=lwd)
points(ts, y[,2], col=col2, pch=20)

# Close device
dev.off()


