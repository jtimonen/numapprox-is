# Run all TMDD experiments

# Sampling
t_start <- Sys.time()
source("sample_bdf.R")
t1 <- Sys.time()

# PSIS
source("reliability.R")
t2 <- Sys.time()

times <- c(t1 - t_start, t2 - t1)
minutes <- round(as.numeric(times / 60), 3)
print(times)
print(minutes)
