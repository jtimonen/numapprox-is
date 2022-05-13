# Helper function
create_gpd_txt <- function(pfit) {
  k_val <- round(pfit$k, 3)
  s_val <- round(pfit$sigma, 3)
  paste0("k = ", k_val, ", ", expression(sigma), " = ", s_val)
}

# Helper function
create_ylab <- function(gq_low, gq_high) {
  t1 <- paste0("low: ", gq_low$solver$abs_tol, ", ", sep = "")
  t2 <- paste0("high: ", gq_high$solver$abs_tol, sep = "")
  paste(t1, t2)
}

# Creates PSIS figure for one reliability run
create_psis_figure <- function(out) {
  rel <- out$reliab

  # Compute importance ratios
  gq_low <- rel$base
  idx_high <- length(rel$files)
  gq_high <- readRDS(file = rel$files[idx_high])
  log_r <- odemodeling::log_ratios(gq_low, gq_high)
  ratios <- exp(as.vector(log_r))
  df <- data.frame(ratios)

  # Determine tail samples
  S <- length(ratios)
  r_eff <- odemodeling::psis_relative_eff(gq_low, gq_high)
  n_tail <- loo:::n_pareto(r_eff, S)
  tail_inds <- (S + 1 - n_tail):S
  ratios_sorted <- sort(ratios)
  ratios_tail <- ratios_sorted[tail_inds]
  rcut <- ratios_sorted[min(tail_inds) - 1]

  # Plot A
  xlab <- parse(text = paste0("r^{M~','~~M^{'*'}}"))
  cvals <- c("1" = "firebrick", "cutoff" = "steelblue")
  plt_A <- ggplot(df, aes(x = ratios)) +
    theme_bw() +
    xlab(xlab) +
    geom_density(alpha = 1, fill = "gray50") +
    geom_vline(aes(xintercept = 1, color = "1")) +
    geom_vline(aes(xintercept = rcut, color = "cutoff")) +
    scale_color_manual(values = cvals) +
    theme(legend.position = c(0.6, 0.6), legend.title = element_blank())

  # Fit generalized Pareto distribution
  rt_centered <- ratios_tail - rcut
  pfit <- loo::gpdfit(rt_centered)
  pd <- evmix::dgpd(rt_centered, sigmau = pfit$sigma, xi = pfit$k)
  df_gpd <- data.frame(rt_centered, pd)
  df_gpd <- df_gpd[2:n_tail, ]
  colnames(df_gpd) <- c("ratios", "density")
  df_centered <- data.frame(rt_centered)
  colnames(df_centered) <- c("ratios")

  # Plot B
  xlab <- xlab <- parse(text = paste0("r^{M~','~~M^{'*'}} - cutoff"))
  plt_B <- ggplot(df_centered, aes(x = ratios)) +
    theme_bw() +
    xlab(xlab) +
    geom_density(alpha = 1, fill = "gray50") +
    geom_line(data = df_gpd, aes(x = ratios, y = density, color = "GPD")) +
    scale_color_manual(values = c("GPD" = "orange")) +
    theme(legend.position = c(0.7, 0.8), legend.title = element_blank()) +
    ylab(create_ylab(gq_low, gq_high))

  # Annotate
  txt <- create_gpd_txt(pfit)
  xval <- max(ratios_tail) / 2
  yval <- max(pd / 2)
  plt_B <- plt_B + annotate("text",
    x = xval, y = yval, label = txt,
    color = "orange"
  )

  # Combine plots
  plt <- ggarrange(plt_A, plt_B, nrow = 2)

  # Return
  list(
    plt = plt, pfit = pfit, n_tail = n_tail, rcut = rcut,
    gq_low = gq_low, gq_high = gq_high
  )
}

# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
p1 <- create_psis_figure(results$outputs[[1]])
p2 <- create_psis_figure(results$outputs[[2]])
p3 <- create_psis_figure(results$outputs[[3]])

plt <- ggarrange(p1$plt, p2$plt, p3$plt, labels = "auto", nrow = 1)
ggsave(plt, file = "figures/tmdd_psis.pdf", width = 9.4, height = 4.3)
