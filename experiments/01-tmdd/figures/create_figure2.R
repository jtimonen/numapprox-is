# Requirements
source("../R/utils.R")
source("../R/models.R")
source("../R/data.R")
source("../R/functions.R")
library(odemodeling)
library(posterior)
library(scales)

# Load results
res_dir <- c("results_bdf")
fn <- paste0("reliability", ".rds")
fp <- file.path(res_dir, fn)
results <- readRDS(file = fp)

# Helper function
time_df <- function(result, ylog) {
  fits <- result$res$fits
  reliab <- result$reliab
  idx <- result$idx
  create_time_comparison_df(fits, reliab, idx, ylog)
}

# Create plot
ylog <- FALSE
df <- NULL
lab0 <- expression(time[MCMC]^{
  BDF(tol)
})
labs <- list()
labs[[1]] <- lab0
for (j in 1:4) {
  out <- results$outputs[[j]]
  tol_bdf <- out$confs[out$idx]
  df_j <- time_df(out, ylog)
  df_j$logtol <- log10(1 / df_j$inv_tol)
  if (j > 1) {
    df_j <- df_j[which(df_j$procedure != "high"), ]
  }
  df_j$procedure <- as.character(df_j$procedure)
  df_j$procedure[which(df_j$procedure != "high")] <- paste0("low", j)
  df <- rbind(df, df_j)
  str <- paste0("time[MCMC]^{BDF(", tol_bdf, ")} + time[PSIS]^{BDF(tol)}")
  labs[[j + 1]] <- parse(text = str)
}
df$procedure <- as.factor(df$procedure)

# Plot
n_yticks <- 8
cols <- c("#010101", "#ca0020", "#f4a582", "#92c5de", "#0571b0")
aesth <- aes(x = logtol, y = time, group = procedure, color = procedure)
plt_A <- ggplot(df, aesth) +
  ylim(0, 10^4) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = cols,
    labels = labs
  ) +
  theme(legend.position = c(0.2, 0.75), legend.title = element_blank()) +
  scale_x_reverse(breaks = unique(round(df$logtol))) +
  xlab("log10(tol)") +
  ylab("time (s)") +
  scale_y_continuous(
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10^x,
      n = n_yticks
    ),
    labels = trans_format("log10", math_format(10^.x))
  )

# Combine
plt <- plt_A
ggsave(plt, filename = "figures/tmdd_figure2.pdf", width = 5.6, height = 3.9)
