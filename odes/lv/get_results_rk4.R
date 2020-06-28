#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

N_CHAINS <- 30
K         <- 9
D         <- 10
STEP_SIZE <- array(0, c(D, K))
PARETO_K  <- array(0, c(D, K))
RUNTIME   <- array(0, c(D, K, N_CHAINS))

for(data_idx in 1:D){
  
  # Read results for this dataset
  for(i in 1:K){
    fn <- paste0('res/rk4/rk4_dat_', data_idx, '_step_', i ,'.rds')
    cat(paste0('Reading file  ', fn, '... '))
    tryCatch({
      res <- readRDS(fn)
      STEP_SIZE[data_idx, i] <- res$step_size
      PARETO_K[data_idx, i]  <- res$pareto_k
      RUNTIME[data_idx, i, ] <- res$runtimes
      cat('success!\n')
    }, error = function(e) {
      STEP_SIZE[data_idx, i] <- NA
      PARETO_K[data_idx, i]  <- NA
      RUNTIME[data_idx, i, ] <- rep(NA, N_CHAINS)
      cat('failed.\n')
    })
  }
}

# Pareto K plot
df1 <- data.frame(as.factor(as.vector(STEP_SIZE)), as.vector(PARETO_K))
colnames(df1) <- c('step_size', 'pareto_k')
p1 <- ggplot(df1, aes(x=step_size, group=step_size, y=pareto_k)) + geom_hline(yintercept=0.5, lty=2, col='firebrick') + 
  geom_boxplot(col='gray20', fill='firebrick2') + 
  theme_bw() + scale_x_discrete(labels=as.character(STEP_SIZE[1,]))
p1 <- p1 + xlab('Step size') + ylab('Pareto k')
p1 <- p1 + scale_y_continuous(breaks=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0))

# Runtime plot
SS <- array(0, dim(RUNTIME))
for(s in 1:dim(SS)[3]){
  SS[,,s] <- STEP_SIZE
}
df2 <- data.frame(as.factor(as.vector(SS)), as.vector(RUNTIME))
colnames(df2) <- c('step_size', 'run_time')
p2 <- ggplot(df2, aes(x=step_size, group=step_size, y=run_time)) +
  geom_boxplot(col='gray20', fill='firebrick2') + 
  theme_bw() + scale_x_discrete(labels=as.character(STEP_SIZE[1,]))
p2 <- p2 + xlab('Step size') + ylab('Run time per chain (s)')
p2 <- p2 + scale_y_continuous(breaks=seq(0,12,by=2))

# Combine plots
plt <- ggarrange(p1, p2, labels='auto')
