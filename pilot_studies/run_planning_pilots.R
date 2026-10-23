# Run or reuse the precision and RF pilots, then update the planning summary.

run_planning_pilots <- function() {
  source("R/workflow_helpers.R")
  sys.source("validate.R", envir = new.env(parent = globalenv()))
  scripts <- c("run/02_run_precision_pilot.R",
               "run/04_run_rf_tree_sensitivity.R",
               "run/10_summarize_planning.R",
               "run/11_export_results.R")
  for (script in scripts) {
    message("Running ", script)
    sys.source(script, envir = new.env(parent = globalenv()))
  }
  return(invisible(TRUE))
}
run_planning_pilots()
