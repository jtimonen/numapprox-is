# Original plot
tol1 <- solver1$abs_tol
tol2 <- solver2$abs_tol
tols <- as.factor(rep(c(tol1, tol2), each = L))
df <- data.frame(as.vector(times), as.vector(ad_calls), tols)
colnames(df) <- c("time", "ad_calls", "tol")
plt <- ggplot(df, aes(x = ad_calls, y = time, group = tol, color = tol)) +
  geom_smooth(
    method = "lm", aes(x = ad_calls, y = time), inherit.aes = FALSE,
    color = "black", se = FALSE, lwd = 0.5
  ) +
  theme_bw() +
  xlab("Number of ODE RHS calls with autodiff") +
  ylab("time (s)") +
  geom_point() +
  scale_color_brewer(type = "qual", palette = 6)

ggsave(plt, filename = "figures/tmdd_profiling.pdf", width = 5.7, height = 4)

# More detailed plot
plt <- ggplot(data = DF, aes(x = tols, y = vals, group = tols, color = tols)) +
  geom_boxplot() +
  facet_wrap(. ~ names, scales = "free_y") +
  ylab("value") +
  xlab(expression(epsilon)) +
  theme_bw() +
  theme(legend.title = element_blank())


ggsave(plt,
  filename = "figures/tmdd_profiling_more.pdf",
  width = 7.5, height = 5.5
)
