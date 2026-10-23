# Reproduce all planning pilots and checks. This does not run the main study.

run_all <- function() {
  source("R/workflow_helpers.R")
  source("R/reproducibility.R")
  check_reference_environment(packages = c("glmnet", "ranger", "batchtools", "snow",
     "dplyr", "tibble", "ggplot2", "Matrix", "xtable"))
  sys.source("validate.R", envir = new.env(parent = globalenv()))
  scripts <- sort(list.files("run", pattern = "^0[1-9]_.*[.]R$", full.names = TRUE))
  scripts <- c(scripts, "run/10_summarize_planning.R", "run/11_export_results.R")
  for (script in scripts) {
    message("Running ", script)
    sys.source(script, envir = new.env(parent = globalenv()))
  }
  return(invisible(TRUE))
}
run_all()
