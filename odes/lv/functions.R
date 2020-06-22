# Create additional Stan data related to some ODE solver methods
add_interpolation_data <- function(data_list, h){
  t0 <- data_list$t0
  ts <- data_list$ts
  R  <- compute_R(t0, ts, h)
  A  <- compute_A(t0, ts, h, R)
  data_list$STEP_SIZE <- h
  data_list$INTERP_R  <- R
  data_list$INTERP_A  <- pmax(A, 0) # rounding errors can make A slightly negative though should be 0
  return(data_list)
}

# Compute integers r_1, ..., r_n
compute_R <- function(t0, ts, h){
  n <- length(ts)
  R <- rep(0, n)
  for(i in 1:n){
    r <- 0
    while(t0 + r*h < ts[i]){
      r <- r + 1
    }
    R[i] <- r - 1
  }
  return(R)
}

# Compute multipliers a_1, ..., a_n
compute_A <- function(t0, ts, h, R){
  n <- length(ts)
  A <- rep(0, n)
  for(i in 1:n){
    D_i <- ts[i] - (t0 + R[i]*h)
    A[i] <- (h - D_i)/h
  }
  return(A)
}

# Helper function
get_samples <- function(stan_fit, param){
  samples <- as.vector(rstan::extract(stan_fit, pars=param)[[param]])
  return(samples)
}
