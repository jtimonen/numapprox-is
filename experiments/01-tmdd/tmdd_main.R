# Run all TMDD experiments
# - requires that simulated_data.rds exists
# - result .rds files will be created in 'results_bdf' directory
# - result figures will be created in 'figures' directory

# MCMC sampling, takes long
source("tmdd_mcmc.R")

# PSIS, is quick
source("tmdd_reliability.R")

# Create result figures
source("figures/create_figure1.R")
