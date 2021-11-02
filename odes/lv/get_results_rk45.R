#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

N_CHAINS <- 100
K <- 7
D <- 10
TOLERANCE <- array(0, c(D, K))
PARETO_K <- array(0, c(D, K))
RUNTIME <- array(0, c(D, K, N_CHAINS))


for (data_idx in 1:D) {

  # Read results for this dataset
  for (i in 1:K) {
    fn <- paste0("res/rk45/rk45_dat_", data_idx, "_tol_", i, ".rds")
    cat(paste0("Reading file  ", fn, "... "))
    tryCatch(
      {
        res <- readRDS(fn)
        TOLERANCE[data_idx, i] <- res$tol
        PARETO_K[data_idx, i] <- res$pareto_k
        RUNTIME[data_idx, i, ] <- res$runtimes
        cat("success!\n")
      },
      error = function(e) {
        TOLERANCE[data_idx, i] <- NA
        PARETO_K[data_idx, i] <- NA
        RUNTIME[data_idx, i, ] <- rep(NA, N_CHAINS)
        cat("failed.\n")
      }
    )
  }
}
TOLERANCE <- log10(TOLERANCE)

# Pareto K plot
df1 <- data.frame(as.factor(as.vector(TOLERANCE)), as.vector(PARETO_K))
colnames(df1) <- c("tol", "pareto_k")
p1 <- ggplot(df1, aes(x = tol, group = tol, y = pareto_k)) +
  geom_hline(yintercept = 0.5, lty = 2, col = "firebrick") +
  geom_boxplot(col = "gray20", fill = "firebrick2") +
  theme_bw() +
  scale_x_discrete(labels = as.character(TOLERANCE[1, ]))
p1 <- p1 + xlab(expression(log[10](tol))) + ylab("Pareto k")

# Runtime plot
TOL <- array(0, dim(RUNTIME))
for (s in 1:dim(TOL)[3]) {
  TOL[, , s] <- TOLERANCE
}
df2 <- data.frame(as.factor(as.vector(TOL)), as.vector(RUNTIME))
colnames(df2) <- c("tol", "run_time")
p2 <- ggplot(df2, aes(x = tol, group = tol, y = run_time)) +
  geom_boxplot(col = "gray20", fill = "firebrick2") +
  theme_bw() +
  scale_x_discrete(labels = as.character(TOLERANCE[1, ]))
p2 <- p2 + xlab(expression(log[10](tol))) + ylab("Run time per chain (s)")
p2 <- p2 + scale_y_continuous(breaks = seq(0, 60, by = 10))

# Combine plots
plt <- ggarrange(p1, p2, labels = "auto")
