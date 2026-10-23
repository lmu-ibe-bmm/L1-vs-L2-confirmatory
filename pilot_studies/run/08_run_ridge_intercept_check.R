# Reproduce the twelve-setting Ridge follow-up.

source("R/workflow_helpers.R")
source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/ridge_intercept_check.R")
ridge_intercept_run <- run_ridge_intercept_check(
  R = 50L, n_test = 5000L, n_cores = project_config$n_cores,
  seed_base = 1700000000L,
  output_dir = "outputs/ridge_intercept_check",
  registry_dir = "outputs/ridge_intercept_check/registry_R50",
  resume = TRUE, overwrite_registry = FALSE,
  resubmit_errors = project_config$resubmit_errors, wait = TRUE, save_outputs = TRUE
)
