# Reproduce the EN follow-up with alpha = 0.95 and 1/3.

source("R/workflow_helpers.R")
source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/elastic_net_intercept_check.R")
elastic_net_intercept_run <- run_elastic_net_intercept_check(
  R = 50L, n_test = 5000L, n_cores = project_config$n_cores,
  seed_base = 1600000000L,
  output_dir = "outputs/elastic_net_intercept_check",
  registry_dir = "outputs/elastic_net_intercept_check/registry_R50_two_alphas",
  resume = TRUE, overwrite_registry = FALSE,
  resubmit_errors = project_config$resubmit_errors, wait = TRUE, save_outputs = TRUE
)
