#!/usr/bin/env Rscript
library(pracma) # for errorbar plot

for(data_idx in 1:10){

  # Read results for this dataset
  K <- 8
  TOL <- rep(0, K)
  PK  <- rep(0, K)
  TIM <- array(0, c(K, 100))
  for(i in 1:K){
    fn <- paste0('res/rk45_dat_', data_idx, '_tol_', i ,'.rds')
    cat(paste0('Reading file  ', fn, '... '))
    tryCatch({
        res <- readRDS(fn)
        cat('success!\n')
        TOL[i] <- res$tol
        PK[i] <- res$pareto_k
        TIM[i, ] <- res$runtimes
    }, error = function(e) {
        cat('failed.\n')
        TOL[i] <- NA
        PK[i] <- NA
        TIM[i] <- NA
    })
  }

  # Create figure
  fig_name <- paste0('figs/rk45_dat', data_idx, '.pdf')
  pdf(fig_name, width = 10, height = 4.6)
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
  plot(x, m, ylab='Runtime per chain (s)', xlab='log10(tol)', 
     ylim=c(0, 30), pch=16, bty="n")
  grid()
  lines(x, m)
  pracma::errorbar(x, m, yerr=s, add=TRUE, bar.col='firebrick',
         bar.len = 0.03)
  points(x, m, pch=16)

  # Close graphics device
  dev.off()
}


