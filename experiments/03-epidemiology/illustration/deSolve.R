library(deSolve)
library(ggplot2)
library(tidyr)
library(ggpubr)

times <- seq(0, 130, by = 0.25)
N_pop <- 1000000

# Basic SEIR system
SEIR <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dS <- -beta * S * I / N_pop
    dE <- beta * S * I / N_pop - a * E
    dI <- a * E - gamma * I
    dR <- gamma * I
    res <- c(dS, dE, dI, dR)
    list(res)
  })
}

forcing <- function(t, eta, xi, t1, nu) {
  eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu)))
}

# SEIR system with control measures
SEIRforcing <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    force <- forcing(t, eta, xi, t1, nu)
    beta_star <- beta * force
    dS <- -beta_star * S * I / N_pop
    dE <- beta_star * S * I / N_pop - a * E
    dI <- a * E - gamma * I
    dR <- gamma * I
    res <- c(dS, dE, dI, dR)
    list(res, force = force)
  })
}

parms <- c(beta = 2, a = 0.2, gamma = 0.7)
parms_f <- c(
  beta = 2, a = 0.2, gamma = 0.7, eta = 0.2,
  nu = 7.6, p_reported = 0.03, t1 = 42, xi = 0.5
)

e0 <- 1
i0 <- 1
xstart <- c(S = N_pop - e0 - i0, E = e0, I = i0, R = 0)

out <- ode(y = xstart, times = times, func = SEIR, parms)
df <- data.frame(out)
df_long <- tidyr::pivot_longer(df, cols = c("S", "E", "I", "R"))
colnames(df_long) <- c("Day", "Group", "Value")

plt_A <- ggplot(df_long, aes(x = Day, y = Value, group = Group, color = Group)) +
  geom_line() +
  ylab("Number of individuals")


out <- ode(y = xstart, times = times, func = SEIRforcing, parms_f)
df <- data.frame(out)
df_seir <- df[, 1:5]
df_long <- tidyr::pivot_longer(df_seir, cols = c("S", "E", "I", "R"))
colnames(df_long) <- c("Day", "Group", "Value")

plt_B <- ggplot(df_long, aes(x = Day, y = Value, group = Group, color = Group)) +
  geom_line() +
  ylab("Number of individuals")

plt <- ggpubr::ggarrange(plt_A, plt_B,
  nrow = 1,
  labels =
    c(
      "no control measures",
      "with control measures"
    )
)

a <- ggplot(df, aes(x = time, y = force)) +
  geom_line(lwd = 1) +
  xlab("t") +
  ylab("f(t)") +
  ylim(0, 1) +
  geom_vline(
    xintercept = 42,
    col = "firebrick", lty = 2
  ) +
  geom_vline(xintercept = 42 + 7.6, col = "steelblue", lty = 2) +
  geom_hline(yintercept = 0.2, color = "orange", lty = 2)
