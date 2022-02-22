# Setup
par_dir <- file.path(res_dir, "rk4")
save_dir <- file.path(par_dir, "mcmc")
odemodeling:::create_dir_if_not_exist(par_dir)
odemodeling:::create_dir_if_not_exist(save_dir)

# Run sampling
num_steps <- 2:20
solvers <- rk4_list(num_steps = num_steps)
fits <- model$sample_manyconf(
  t0 = dat$t0, t = dat$t, data = add_data, init = init,
  solvers = solvers, step_size = 0.1, savedir = save_dir,
  iter_warmup = ITER, iter_sampling = ITER, chains = CHAINS
)

# Save results
results <- list(
  fits = fits,
  session_info = sessionInfo()
)
fp <- file.path(par_dir, "mcmc.rds")
saveRDS(results, file = fp)
