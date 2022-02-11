
# Plot incidence
df_dat <- data.frame(dat$cases, dat$ts)
colnames(df_dat) <- c("cases", "t")
df <- get_incidence_quantiles_df(fit)
map <- aes(x = t, y = cases)
plt <- ggplot(df, aes(x = t, y = median, ymin = lower, ymax = upper)) +
  geom_line(color = "firebrick") +
  geom_ribbon(fill = "firebrick2", alpha = 0.6) +
  geom_point(data = df_dat, mapping = map, inherit.aes = FALSE) +
  ylab("Number of reported cases")
