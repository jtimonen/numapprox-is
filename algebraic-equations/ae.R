
library(rstan)
n   <- 300
sigma <- 0.5
x1  <- rnorm(n, 1, sigma)
x2  <- rnorm(n, -1, sigma)
x   <- cbind(x1, x2)
plot(x, xlim=c(-3,3), ylim=c(-3,3))
grid()

rtol <- 1
max_iter <- 1000
dat <- list(x = x, n = n, sigma = sigma, rtol=rtol, max_iter=max_iter)
sm  <- stan_model(file='ae.stan')
fit <- sampling(sm, data = dat, chains = 4, cores = 1, 
                control=list(adapt_delta=0.95), iter = 2000)
print(fit)

