# Run all Lotka-Volterra experiments

# Sampling
t_start <- Sys.time()
source("sample_midpoint.R")
t1 <- Sys.time()
source("sample_rk4.R")
t2 <- Sys.time()
source("sample_rk45.R")
t3 <- Sys.time()

# PSIS
source("reliability.R")
t4 <- Sys.time()

times <- c(t1 - t_start, t2 - t1, t3 - t2, t4 - t3)
minutes <- round(as.numeric(times / 60), 3)
print(times)
print(minutes)
