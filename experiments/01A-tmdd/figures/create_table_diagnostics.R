# Diagnose sampler (note: only 2 significant digits)
diagnose_sampler <- function(fit) {
  fn <- paste0(tempfile(), ".csv")
  cat(fit$cmdstan_summary$stdout, file = fn)
  r <- read.csv(file = fn, sep = "\t")
  r <- as.list(r[6:(6 + 4), ])
  splitter <- function(x) strsplit(x, " +")
  r <- sapply(r, splitter)
  unlink(fn)
  r <- matrix(unlist(r), 5, 10, byrow = T)[, 1:2]
  nams <- r[, 1]
  mean <- r[, 2]
  df <- data.frame(mean)
  rownames(df) <- nams
  df
}

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
fits <- results$outputs[[1]]$res$fits # mcmc results

L <- length(fits$files)
mcmc_diags <- NULL
hmc_diags <- NULL
for (j in seq_len(L)) {
  fit <- readRDS(fits$files[j])
  smr <- fit$summary()[1:9, ] # 8 params + lp
  atol <- fit$solver$abs_tol
  time <- fit$time()$total
  max_rhat <- max(smr$rhat)
  min_eb <- min(smr$ess_bulk)
  min_et <- min(smr$ess_tail)
  a <- diagnose_sampler(fit)
  mcmc_diags <- rbind(mcmc_diags, c(atol, time, max_rhat, min_eb, min_et))
  hmc_diags <- rbind(hmc_diags, c(atol, a$mean))
}
mcmc_diags <- data.frame(mcmc_diags)
colnames(mcmc_diags) <- c(
  "epsilon", "runtime", "max_rhat", "min_ess_bulk", "min_ess_tail"
)
colnames(hmc_diags) <- c("epsilon", row.names(a))

#  Create tables
library(xtable)
t1 <- xtable(mcmc_diags,
  display = c("d", "e", "f", "f", "f", "f"), digits = 3,
  label = "t: mcmc_diagnostics_tmdd",
  caption = "MCMC diagnostics in the TMDD experiment."
)
t2 <- xtable(hmc_diags,
  digits = 2,
  label = "t: hmc_diagnostics_tmdd",
  caption = "HMC NUTS diagnostics in the TMDD experiment."
)
