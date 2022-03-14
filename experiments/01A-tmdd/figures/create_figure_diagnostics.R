# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
fits <- results$outputs[[1]]$res$fits # mcmc results

L <- length(fits$files)
diags <- NULL
for (j in seq_len(L)) {
  fit <- readRDS(fits$files[j])
  smr <- fit$summary()[1:9, ] # 8 params + lp
  atol <- fit$solver$abs_tol
  time <- fit$time()$total
  max_rhat <- max(smr$rhat)
  min_eb <- min(smr$ess_bulk)
  min_et <- min(smr$ess_tail)
  diags <- rbind(diags, c(atol, time, max_rhat, min_eb, min_et))
}
diags <- data.frame(diags)
colnames(diags) <- c("tol", "runtime", "max_rhat", "min_ess_bulk", "min_ess_tail")

# TODO: gradient diagnose
