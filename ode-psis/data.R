# Lynx-hare data ----------------------------------------------------------

# Following https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html
load_data_lynxhare <- function() {
  df <- raw_data_lynxhare()
  N <- length(df$year) - 1
  ts <- 1:N
  y_init <- c(df$hare[1], df$lynx[1])
  y <- as.matrix(df[2:(N + 1), 2:3])
  y <- cbind(y[, 2], y[, 1]) # hare, lynx order
  colnames(y) <- c("hare", "lynx")

  # Return
  list(t0 = df$year[1], t = df$year[2:(N + 1)], y_obs_init = y_init, y_obs = y)
}

# Data from http://www.math.tamu.edu/~phoward/m442/modbasics.pdf
# Downloaded 15 October 2017, 4:59 PM EDT
raw_data_lynxhare <- function() {

  # Year, Lynx, Hare
  x <- c(
    1900, 4.0, 30.0,
    1901, 6.1, 47.2,
    1902, 9.8, 70.2,
    1903, 35.2, 77.4,
    1904, 59.4, 36.3,
    1905, 41.7, 20.6,
    1906, 19.0, 18.1,
    1907, 13.0, 21.4,
    1908, 8.3, 22.0,
    1909, 9.1, 25.4,
    1910, 7.4, 27.1,
    1911, 8.0, 40.3,
    1912, 12.3, 57.0,
    1913, 19.5, 76.6,
    1914, 45.7, 52.3,
    1915, 51.1, 19.5,
    1916, 29.7, 11.2,
    1917, 15.8, 7.6,
    1918, 9.7, 14.6,
    1919, 10.1, 16.2,
    1920, 8.6, 24.7
  )
  X <- matrix(x, ncol = 3, byrow = TRUE)
  df <- data.frame(X)
  colnames(df) <- c("year", "lynx", "hare")
  return(df)
}
