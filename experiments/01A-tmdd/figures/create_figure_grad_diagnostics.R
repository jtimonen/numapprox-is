# Euclidean norm of vector
euc_norm <- function(x) {
  sqrt(sum(x * x))
}

subsample_draws <- function(fit, draw_inds) {
  draws <- posterior::merge_chains(fit$draws())
  pnames <- c("k_on", "k_off", "k_in", "k_out", "k_eL", "k_eP", "sigma")
  posterior:::subset_draws(
    draws,
    iteration = draw_inds,
    variable = pnames
  )
}

create_param_list <- function(fit, draw_ind) {
  draws <- posterior::merge_chains(fit$draws())
  pnames <- c("k_on", "k_off", "k_in", "k_out", "k_eL", "k_eP", "sigma")
  out <- list()
  j <- 0
  for (var in pnames) {
    j <- j + 1
    d <- posterior:::subset_draws(
      draws,
      iteration = draw_ind,
      variable = var
    )
    out[[j]] <- as.numeric(d)
  }
  names(out) <- pnames
  return(out)
}

# Diagnose gradient
diagnose_grad <- function(model, fit1, fit, draw_ind, ...) {
  params <- create_param_list(fit1, draw_ind)
  params <- rep(list(params), 1) # same for all chains
  model$diagnose(
    t0 = fit1$t0, t = fit1$t, data = fit1$data,
    solver = fit$solver, init = params, ...
  )
}

# Compute error in gradient
compute_grad_error <- function(model, fit1, fit, fn_epsilon, draw_ind) {
  g <- diagnose_grad(model, fit1, fit,
    draw_ind = draw_ind,
    epsilon = fn_epsilon
  )
  lp <- g$lp
  grad_fn <- g$gradients$finite_diff
  grad_error <- g$gradients$error
  en_err <- euc_norm(grad_error)
}

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
fits <- results$outputs[[1]]$res$fits # mcmc results
fit1 <- readRDS(fits$files[1])
n_draws <- 100
draw_inds <- sample.int(fit1$ndraws(), n_draws)

fn_epsilon <- 1e-6 # for finite diff
FN_EPS <- c(1e-3, 1e-6, 1e-9)
plots <- list()
for (plt_idx in 1:3) {
  fn_epsilon <- FN_EPS[plt_idx]

  L <- length(fits$files)
  diags <- NULL
  log_errs <- matrix(0, L, n_draws)
  tols <- matrix(0, L, n_draws)
  for (j in seq_len(L)) {
    fit <- readRDS(fits$files[j])
    atol <- fit$solver$abs_tol
    k_idx <- 0
    for (k in draw_inds) {
      k_idx <- k_idx + 1
      err <- compute_grad_error(model, fit1, fit, fn_epsilon, k)
      log_errs[j, k_idx] <- log(err)
      tols[j, k_idx] <- atol
      cat("*")
    }
    cat(" j =", j, "\n")
    diags <- rbind(diags, c(atol, log(err)))
  }
  diags <- data.frame(diags)
  colnames(diags) <- c("tol", "log_err_norm_grad")

  # Create data frame
  tol <- as.numeric(tols)
  df <- data.frame(as.numeric(log_errs), as.factor(tol))
  colnames(df) <- c("grad_error", "tol")
  plt <- ggplot(df, aes(x = tol, group = tol, y = grad_error)) +
    geom_boxplot(fill = "firebrick3") +
    theme_bw() +
    xlab(expression(epsilon)) +
    ggtitle(paste0("delta = ", fn_epsilon))
  plots[[plt_idx]] <- plt
}

turner <- function(x) {
  x + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    ylab("log(grad_error_norm)")
}
plots2 <- lapply(plots, turner)
plt <- ggarrange(plotlist = plots2, labels = "auto", nrow = 1)
ggsave(plt, file = "figures/tmdd_grad_error.pdf", width = 10, height = 3.45)
