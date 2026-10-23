# Run the shared-core, planning and diagnostic tests.
# Unexpected warnings stop the tests. Version differences remain warnings.

run_all_tests <- function() {
  old_warn <- getOption("warn")
  options(warn = 2)
  on.exit(options(warn = old_warn), add = TRUE)
  files <- c("test_reproducibility.R", "test_core.R", "test_precision_pilot.R",
             "test_runtime_pilot.R", "test_testset_size_pilot.R",
             "test_rf_tree_sensitivity.R", "test_signal_sanity.R",
             "test_lasso_intercept_check.R", "test_elastic_net_intercept_check.R",
             "test_ridge_intercept_check.R", "test_workflow.R")
  for (file in files) {
    cat("\nRunning:", file, "\n")
    withCallingHandlers(
      source(file.path("tests", file), local = new.env(parent = globalenv())),
      reference_version_warning = function(condition) {
        message(conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    cat("OK:", file, "\n")
  }
  cat("\nALL PROJECT TESTS PASSED.\n")
  return(invisible(TRUE))
}
run_all_tests()
