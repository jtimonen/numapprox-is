# Setup
par_dir <- file.path(res_dir, "rk45")
save_dir <- file.path(par_dir, "mcmc")
odemodeling:::create_dir_if_not_exist(par_dir)
odemodeling:::create_dir_if_not_exist(save_dir)

# Run sampling
tols <- c(0.0075, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 10^(-5:-12))
solvers <- rk45_list(tols = tols, max_num_steps = 1e5)
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
