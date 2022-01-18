# Colorize string
colorize_str <- function(x, col = "\u001b[34m") {
  paste0(col, x, "\u001b[0m")
}

# Create numeric index for the experiment based on command line args
setup_experiment_index <- function(args) {
  if (interactive()) {
    idx <- 0
  } else {
    idx <- args[1]
  }
  idx <- as.numeric(idx)
  msg <- paste("\n-------------- idx = ", idx, " ----------------\n", sep = "")
  cat(colorize_str(msg))
  return(idx)
}

# Setup paths and result file names
setup_experiment_paths <- function(idx, res_dir = "res") {
  odemodeling:::create_dir_if_not_exist(res_dir)
  res_file <- file.path(res_dir, paste0("res_", idx, ".rds"))
  message("results will be saved to: ", res_file, sep = "")
  list(
    res_file = res_file,
    res_dir = res_dir
  )
}
