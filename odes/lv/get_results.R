library(pracma) # for errorbar plot

K <- 8
TOL <- rep(0, K)
PK  <- rep(0, K)
TIM <- array(0, c(K, 100))
for(i in 1:K){
    fn <- paste0('res/res_rk45_', i ,'.rds')
    cat(paste0('Reading file  ', fn, '\n'))
    res <- readRDS(fn)
    TOL[i] <- res$tol
    PK[i] <- res$pareto_k
    TIM[i, ] <- res$runtimes
}

par(mfrow=c(1,2))
# Pareto K plot
x <- log10(TOL)
plot(x, PK, ylab='Pareto k', xlab='log10(tol)', pch=16, bty="n")
grid()
lines(x, PK)
points(x, PK, pch=16)

# Runtime plot
m <- rowMeans(TIM)
s <- apply(TIM, 1, stats::sd)
plot(x, m, ylab='Runtime (s)', xlab='log10(tol)', 
     ylim=c(0, 30), pch=16, bty="n")
grid()
lines(x, m)
pracma::errorbar(x, m, yerr=s, add=TRUE, bar.col='firebrick',
         bar.len = 0.03)
points(x, m, pch=16)
