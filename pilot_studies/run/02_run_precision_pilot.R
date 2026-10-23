# Re-estimate the precision requirement over all 192 settings.

source("R/workflow_helpers.R")
source("R/precision_pilot.R")

output_dir <- "outputs/precision_pilot_R500"
source_files <- c("R/pilot_core.R", "R/reproducibility.R", "R/precision_pilot.R")
specification <- list(R = 500L, seed_base = project_config$precision_seed_base,
                      n_settings = 192L, n_test = 5000L,
                      multipliers = c(1, 1.1, 1.2, 1.3), thetas = c(-0.02, 0, 0.02))
precision_pilot <- load_completed_run(output_dir, specification, source_files)

if (is.null(precision_pilot)) {
  message("Precision: 192 settings x 500 repetitions = 96,000 tasks.")
  precision_pilot <- run_precision_pilot(
    R_sd = project_config$precision_R,
    output_dir = output_dir,
    backend = "batchtools",
    n_cores = project_config$n_cores,
    seed_base = project_config$precision_seed_base,
    registry_dir = file.path(output_dir, "registry_R500"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = project_config$resubmit_errors,
    wait = TRUE,
    save_outputs = TRUE,
    sd_multiplier_grid = specification$multipliers,
    theta_assumed_grid = specification$thetas,
    source_files = normalizePath("R/pilot_core.R")
  )
  check_precision_run(precision_pilot)
  precision_pilot$pilot_repetitions <- 500L
  precision_pilot <- save_run_session(precision_pilot, output_dir)
  mark_run_completed(precision_pilot, output_dir, specification, source_files)
}
check_precision_run(precision_pilot)
print(precision_pilot$global_summary)
print(precision_pilot$sensitivity_summary)
