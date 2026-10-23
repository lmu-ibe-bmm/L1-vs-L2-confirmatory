# Small signal and RF checks.

cat("\n=== Signal sanity integration tests ===\n")

if (!file.exists("R/pilot_core.R") || !file.exists("R/precision_pilot.R") || !file.exists("R/signal_sanity.R")) {
  stop("Run tests from the project root after unpacking the signal-sanity package.", call. = FALSE)
}

source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/signal_sanity.R")

signal_sanity_assert_core()
stopifnot(requireNamespace("ranger", quietly = TRUE))
stopifnot(requireNamespace("batchtools", quietly = TRUE))

beta_test <- c(1, 2, rep(0, 8))
v0 <- signal_sanity_analytic_signal_variance(beta_test, rho = 0, b = 10L)
v5 <- signal_sanity_analytic_signal_variance(beta_test, rho = 0.5, b = 10L)
stopifnot(abs(v0 - 5) < 1e-12)
stopifnot(abs(v5 - (0.5 * 5 + 0.5 * 9)) < 1e-12)

design_test <- signal_sanity_prepare_design(n_test = 500L)
settings_test <- signal_sanity_make_settings(design_test, expected_n_settings = 192L)
tasks_test <- signal_sanity_make_tasks(settings_test, R = 50L, seed_base = 1300000000L)
stopifnot(nrow(settings_test) == 192L)
stopifnot(nrow(tasks_test) == 9600L)
stopifnot(!anyDuplicated(tasks_test$seed))

local({
  fault_worker <- signal_sanity_one_task
  fault_environment <- new.env(parent = environment(fault_worker))
  fault_environment$simulate_data <- function(...) stop("Controlled data-generation failure.")
  environment(fault_worker) <- fault_environment
  fault_task <- tasks_test[1L, , drop = FALSE]
  failed_result <- fault_worker(task_row = fault_task, settings = settings_test, design = design_test, rf_num_trees = 5L, analysis_n_test = 500L)
  stopifnot(isFALSE(failed_result$ok[[1L]]),
            grepl("Controlled data-generation failure", failed_result$error_message[1L], fixed = TRUE))
  for (field in c("task_id", "setting_id", "run_id", "seed")) {
    stopifnot(identical(as.integer(failed_result[[field]][1L]),
                        as.integer(fault_task[[field]][1L])))
  }
  fault_task$seed <- NULL
  malformed_result <- fault_worker(task_row = fault_task, settings = settings_test, design = design_test, rf_num_trees = 5L, analysis_n_test = 500L)
  stopifnot(isFALSE(malformed_result$ok[[1L]]),
            is.na(malformed_result$seed[1L]),
            malformed_result$task_id[1L] == fault_task$task_id[1L],
            malformed_result$setting_id[1L] == fault_task$setting_id[1L],
            malformed_result$run_id[1L] == fault_task$run_id[1L])
})
cat("Structured failure identifiers passed.\n")

one_task <- signal_sanity_one_task(
  task_row = tasks_test[1L, , drop = FALSE],
  settings = settings_test,
  design = design_test,
  rf_num_trees = 5L,
  analysis_n_test = 500L
)
stopifnot(nrow(one_task) == 1L)
stopifnot(isTRUE(one_task$ok[[1L]]))
stopifnot(is.finite(one_task$mse_intercept[[1L]]), one_task$mse_intercept[[1L]] > 0)
stopifnot(is.finite(one_task$mse_oracle[[1L]]), one_task$mse_oracle[[1L]] > 0)
stopifnot(is.finite(one_task$mse_rf[[1L]]), one_task$mse_rf[[1L]] > 0)
stopifnot(is.finite(one_task$analytic_snr[[1L]]))
stopifnot(abs(one_task$analytic_snr_rel_error[[1L]]) < 1e-6)
stopifnot(one_task$n_test[[1L]] == 500L)
stopifnot(one_task$n_test_generated[[1L]] >= one_task$n_test[[1L]])
stopifnot(one_task$rf_num_trees[[1L]] == 5L)

cat("End-to-end test setting: ", one_task$setting_id, "\n", sep = "")
cat("Analytic SNR relative error: ", format(one_task$analytic_snr_rel_error, scientific = TRUE), "\n", sep = "")
cat("Intercept MSE: ", signif(one_task$mse_intercept, 5), "\n", sep = "")
cat("Oracle MSE:    ", signif(one_task$mse_oracle, 5), "\n", sep = "")
cat("RF MSE (5-tree API test only): ", signif(one_task$mse_rf, 5), "\n", sep = "")
cat("=== Direct integration test passed ===\n\n")

cat("=== Batchtools/socket smoke test ===\n")
bt_smoke <- signal_sanity_batchtools_smoke_test(
  settings = settings_test,
  design = design_test,
  task_row = tasks_test[1L, , drop = FALSE],
  analysis_n_test = 500L,
  rf_num_trees = 5L,
  project_root = ".",
  n_cores = 1L
)
stopifnot(nrow(bt_smoke) == 1L)
stopifnot(isTRUE(bt_smoke$ok[[1L]]))
cat("Batchtools/socket smoke test passed.\n")
cat("=== ALL integration tests passed ===\n\n")
