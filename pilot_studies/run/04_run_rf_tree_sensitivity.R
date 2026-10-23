# Repeat the tree-count check with corrected, separate seeds.

source("R/workflow_helpers.R")
source("R/rf_tree_sensitivity.R")

output_dir <- "outputs/rf_tree_sensitivity_R100"
source_files <- c("R/pilot_core.R", "R/reproducibility.R",
                  "R/precision_pilot.R", "R/rf_tree_sensitivity.R")
specification <- list(R = 100L, seed_base = project_config$rf_seed_base,
                      seed_policy = "separated", n_settings = 192L,
                      tree_counts = c(200L, 500L, 1000L),
                      multipliers = c(1, 1.1, 1.2, 1.3), thetas = c(-0.02, 0, 0.02))
rf_tree_sensitivity <- load_completed_run(output_dir, specification, source_files)

if (is.null(rf_tree_sensitivity)) {
  message("RF: 192 settings x 100 repetitions = 19,200 tasks.")
  rf_tree_sensitivity <- run_rf_tree_sensitivity(
    R_sensitivity = project_config$rf_R,
    tree_counts = specification$tree_counts,
    setting_ids = NULL,
    output_dir = output_dir,
    backend = "batchtools",
    n_cores = project_config$n_cores,
    seed_base = project_config$rf_seed_base,
    seed_policy = "separated",
    registry_dir = file.path(output_dir, "registry_R100_separated"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = project_config$resubmit_errors,
    wait = TRUE,
    save_outputs = TRUE,
    sd_multiplier_grid = specification$multipliers,
    theta_assumed_grid = specification$thetas,
    source_files = normalizePath("R/pilot_core.R")
  )
  check_rf_run(rf_tree_sensitivity)
  rf_tree_sensitivity <- save_run_session(rf_tree_sensitivity, output_dir)
  mark_run_completed(rf_tree_sensitivity, output_dir, specification, source_files)
}
check_rf_run(rf_tree_sensitivity)
print(rf_tree_sensitivity$tree_count_summary)
print(rf_tree_sensitivity$sensitivity_summary)
