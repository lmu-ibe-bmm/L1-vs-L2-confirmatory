# Read the two job registries without resetting anything.

source("R/workflow_helpers.R")
paths <- c(precision = "outputs/precision_pilot_R500/registry_R500",
           rf = "outputs/rf_tree_sensitivity_R100/registry_R100_separated")
for (name in names(paths)) {
  cat("\n", name, "\n", sep = "")
  if (!file.exists(file.path(paths[[name]], "registry.rds"))) {
    output_dir <- dirname(paths[[name]])
    record_file <- file.path(output_dir, "run_complete.rds")
    if (file.exists(record_file)) {
      record <- readRDS(record_file)
      result <- load_completed_run(output_dir, record$specification,
                                   names(record$source_md5))
      if (is.null(result)) stop("Completed result files are missing: ", output_dir)
      cat("Completed results are available. No job registry is included here.\n")
    } else {
      cat("No job registry or completed results found.\n")
    }
    next
  }
  reg <- batchtools::loadRegistry(paths[[name]], writeable = FALSE, make.default = FALSE)
  print(batchtools::getStatus(reg = reg))
  errors <- batchtools::findErrors(reg = reg)
  if (nrow(errors)) {
    cat("Failed jobs:\n")
    print(errors)
    print(batchtools::getErrorMessages(ids = errors, reg = reg))
  }
}
