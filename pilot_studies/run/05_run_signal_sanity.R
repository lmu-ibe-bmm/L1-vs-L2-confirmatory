# Reproduce the full-grid signal and RF check.

source("R/workflow_helpers.R")
source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/signal_sanity.R")
signal_sanity_run <- run_signal_sanity_full_grid(
  R = 50L, n_test = 5000L, rf_num_trees = 200L,
  n_cores = project_config$n_cores, seed_base = 1300000000L,
  output_dir = "outputs/signal_sanity",
  registry_dir = "outputs/signal_sanity/registry_R50",
  resume = TRUE, overwrite_registry = FALSE,
  resubmit_errors = project_config$resubmit_errors, wait = TRUE, save_outputs = TRUE
)
