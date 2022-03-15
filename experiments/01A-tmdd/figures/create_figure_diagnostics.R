# Euclidean norm of vector
euc_norm <- function(x) {
  sqrt(sum(x * x))
}

# Diagnose gradient
diagnose_grad <- function(model, fit, ...) {
  model$diagnose(
    t0 = fit$t0, t = fit$t, data = fit$data,
    solver = fit$solver, ...
  )
}

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
fits <- results$outputs[[1]]$res$fits # mcmc results

fn_epsilon <- 1e-6 # for finite diff
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
  g <- diagnose_grad(model, fit, init = init, epsilon = fn_epsilon)
  lp <- g$lp
  grad_fn <- g$gradients$finite_diff
  grad_error <- g$gradients$error
  en_err <- euc_norm(grad_error)
  diags <- rbind(diags, c(atol, time, max_rhat, min_eb, min_et, lp, log(en_err)))
}
diags <- data.frame(diags)
colnames(diags) <- c(
  "tol", "runtime", "max_rhat", "min_ess_bulk",
  "min_ess_tail", "lp", "log_err_norm_init_grad"
)


plot(-log10(diags$tol), diags$log_err_norm_init_grad, "o", pch = 16)
grid()
