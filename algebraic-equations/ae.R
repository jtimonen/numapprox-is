
library(rstan)
n   <- 300
sigma <- 0.5
x  <- rnorm(n, 1.5, sigma)
plot(x)
grid()

rel_tol <- 1e-10
max_steps <- 1000
dat <- list(x = x, n = n, sigma = sigma, rel_tol=rel_tol, max_steps=max_steps,
            mu_guess=as.array(c(2.0)))
sm  <- stan_model(file='ae.stan')
fit <- sampling(sm, data = dat, chains = 4, cores = 1,
                control=list(adapt_delta=0.95), iter = 10000)
print(fit)

times <- get_elapsed_time(fit)
t <- mean(rowSums(times))
cat(paste("\nAverage runtime per chain:", t, "s.\n"))
