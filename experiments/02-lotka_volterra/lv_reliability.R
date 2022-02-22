# Define function
lv_reliability <- function(dir, idx) {
  res <- readRDS(file.path(dir, "mcmc.rds"))
  fits <- res$fits
  solver1 <- fits$solvers[[1]]
  if (is(solver1, "AdaptiveOdeSolver")) {
    confs <- get_tol_vec(fits$solvers)
  } else {
    confs <- get_num_steps_vec(fits$solvers)
  }

  # Load fit
  inds_rel <- (1 + idx):length(confs)
  fit <- load_fit(file = fits$files[idx])
  confs_rel <- confs[inds_rel]
  if (solver1$name == "rk45") {
    rel_solvers <- rk45_list(tols = confs_rel)
  } else if (solver1$name == "rk4") {
    rel_solvers <- rk4_list(num_steps = confs_rel)
  } else if (solver1$name == "midpoint") {
    rel_solvers <- midpoint_list(num_steps = confs_rel)
  } else {
    stop("unknown solver")
  }

  # Run reliability check
  reliab <- fit$reliability(
    solvers = rel_solvers, force = TRUE, savedir = dir
  )

  # Return list
  list(
    reliab = reliab,
    confs = confs,
    confs_rel = confs_rel,
    res = res,
    idx = idx
  )
}

# Run
dirs <- c("rk45", "rk4", "midpoint")
par_dirs <- file.path(res_dir, dirs)
inds <- c(3, 1, 1)
for (j in 1:3) {
  par_dir <- par_dirs[j]
  fp <- file.path(par_dir, "reliability.rds")
  out <- lv_reliability(par_dir, inds[j])
  saveRDS(out, file = fp)
}
