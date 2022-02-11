# Lynx-hare data ----------------------------------------------------------

# Following https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html
load_data_lynxhare <- function() {
  df <- odemodeling::lynxhare
  N <- length(df$year) - 1
  ts <- 1:N
  y_init <- c(df$hare[1], df$lynx[1])
  y <- as.matrix(df[2:(N + 1), 2:3])
  y <- cbind(y[, 2], y[, 1]) # hare, lynx order
  colnames(y) <- c("hare", "lynx")

  # Return
  list(
    t0 = df$year[1],
    t = df$year[2:(N + 1)],
    y_obs_init = y_init,
    y_obs = y
  )
}

# Swiss data
load_data_switzerland <- function(parent_dir = ".") {
  fp <- file.path(parent_dir, "switzerland.rds")
  readRDS(file = fp)
}

# Following https://github.com/jriou/covid_adjusted_cfr/
load_data_lombardia <- function(parent_dir = ".") {
  fp <- file.path(parent_dir, "data_lombardia.rds")
  readRDS(file = fp)
}
