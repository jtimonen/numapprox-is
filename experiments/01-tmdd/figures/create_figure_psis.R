# Load results
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)
out <- results$outputs[[1]]
rel <- out$reliab

# Compute importance ratios
gq_low <- rel$base
gq_high <- readRDS(file = rel$files[15])
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
plt_A <- ggplot(df, aes(x = ratios)) + # geom_histogram(bins=80) +
  geom_density(alpha = 0.5, fill = "firebrick") +
  geom_vline(xintercept = rcut, lty = 2, color = "steelblue4") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  xlab(xlab)

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
  geom_density(alpha = 0.5, fill = "gray20") +
  geom_line(data = df_gpd, aes(x = ratios, y = density), color = "firebrick") +
  xlab(xlab)


# Combine
plt <- ggarrange(plt_A, plt_B, labels = "auto", nrow = 1)

ggsave(plt, file = "figures/tmdd_psis.pdf", width = 7.9, height = 2.65)
