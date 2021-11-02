
quintic <- function(x, a) {
  return(-x^5 + x + a)
}

t <- seq(-1.4, 1.4, by = 0.01)
y1 <- quintic(t, 1.0)
y2 <- quintic(t, 1.25)
y3 <- quintic(t, 1.5)
y4 <- quintic(t, 2.0)

plot(t, y1, "l", ylim = c(-1, 3))
lines(t, y2, "l")
lines(t, y3, "l")
lines(t, y4, "l")
grid()
