#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
model_idx <- as.numeric(args[3])
step_size <- as.numeric(args[4])
chains    <- as.numeric(args[1])
iter      <- as.numeric(args[2])

library(rstan)
library(loo)
source("functions.R")

# Read data
data_list_model13 <- readRDS('data13.rds')

# Compile model
if(model_idx == 1){
  file <- 'model13/model13_euler.stan'
}else if(model_idx == 2 ){
  file <- 'model13/model13_emp.stan'
}else if(model_idx == 3 ){
  file <- 'model13/model13_ab2.stan'
}else if(model_idx == 4 ){
  file <- 'model13/model13_rk4.stan'
}else if(model_idx == 5 ){
  file <- 'model13/model13_rk45.stan'
}else if(model_idx == 6 ){
  file <- 'model13/model13_bdf.stan'
}else{
  stop("invalid model_idx!")
}

# Compile model
cat("COMPILING STAN MODEL", file, "\n")
model <- stan_model(file = file)
cat("DONE!\n")

# Additional data for reference method ode_integrate_bdf
data_list_model13$EPS           <- 1.0E-9
data_list_model13$abs_tol_REF_  <- 1.0E-10
data_list_model13$rel_tol_REF_  <- 1.0E-10
data_list_model13$max_iter_REF_ <- 1.0E6

# Additional data needed if using a fixed-step solver
if(model_idx < 5){
   cat("USING STEP SIZE", step_size, "\n")
}else{
  cat("USING AN ADAPTIVE STEP SIZE METHOD\n")
}

data_list_model13 <- add_interpolation_data(data_list_model13, step_size)


# Create file name for saving
ss <- gsub('[.]', '_', as.character(step_size))
save_name <- paste0('res/res-', model_idx, "-", ss, ".rds")
cat("SAVE_NAME IS", save_name, "\n")

cat("RUNNING", chains, "CHAINS FOR", iter, "ITERATIONS\n")

# Run sampling
fit <- sampling(object  = model,
                data    = data_list_model13,
                iter    = iter,
                chains  = chains)

# Extract log posterior values (not Jacobian adjusted)
lh1 <- get_samples(fit, 'log_lik_na')
lh2 <- get_samples(fit, 'log_lik_na_REF_')
pr1 <- get_samples(fit, 'log_prior_na')
pr2 <- get_samples(fit, 'log_prior_na_REF_')
post1 <- lh1 + pr1
post2 <- lh2 + pr2

# PSIS
out <- psis(post1 - post2)
pareto_k <- out$diagnostics$pareto_k
cat('PARETO_K =', pareto_k, '\n')

# Get runtime
runtimes <-get_elapsed_time(fit)
print(runtimes)

# Save result
res <- list(fit=fit, pareto_k=pareto_k, runimes=runtimes)
saveRDS(res, file=save_name)


