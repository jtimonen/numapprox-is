#!/usr/bin/env Rscript
library(rstan)

# Settings
by       <- 0.01
T_max    <- 6
y0       <- c(1, 1)
sigma    <- 1 # has no effect
col1     <- 'firebrick'
col2     <- 'steelblue3'
col3     <- 'orange'

# Compile simulation model
sm <- stan_model(file='stan/lv_simulate.stan')

# Crate figure
par(mfrow=c(1,2))

# Solve 1
th1   <- c(1, 2)
ts <- seq(by, T_max, by=by)
d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th1, sigma=sigma)
f1 <- sampling(sm, data = d1, chains = 1, cores = 1, iter = 1, algorithm="Fixed_param")
y1 <- rstan::extract(f1)$y_hat


# Solve 2
th2   <- c(2, 1)
ts <- seq(by, T_max, by=by)
d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th2, sigma=sigma)
f1 <- sampling(sm, data = d1, chains = 1, cores = 1, iter = 1, algorithm="Fixed_param")
y2 <- rstan::extract(f1)$y_hat
for(j in 1:2){
  
}

# Solve 3
th3   <- c(2, 2)
ts <- seq(by, T_max, by=by)
d1 <- list(T=length(ts), y0=y0, t0=0, ts=ts, theta=th3, sigma=sigma)
f1 <- sampling(sm, data = d1, chains = 1, cores = 1, iter = 1, algorithm="Fixed_param")
y3 <- rstan::extract(f1)$y_hat

# PLot
LWD <- 2
plot(ts, y1[1,,1], 'l', col=col1, ylim = c(0, 5), xlab='t', ylab=expression(y[1](t)), lwd=LWD)
lines(ts, y2[1,,1], 'l', col=col2, lwd=LWD)
lines(ts, y3[1,,1], 'l', col=col3, lwd=LWD)

plot(ts, y1[1,,2], 'l', col=col1, ylim = c(0, 5), xlab='t', ylab=expression(y[2](t)), lwd=LWD)
lines(ts, y2[1,,2], 'l', col=col2, lwd=LWD)
lines(ts, y3[1,,2], 'l', col=col3, lwd=LWD)

legend(2.5, 5, c(expression(psi ~"=(1,2)"), expression(psi ~'=(2,1)'), expression(psi~"=(2,2)")),
               lwd=c(2,2,2), col=c(col1, col2, col3), bty="n", cex=0.8)

