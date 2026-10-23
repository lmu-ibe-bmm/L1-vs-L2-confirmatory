# Small LASSO fits and collected gains.

cat("\n=== LASSO / INTERCEPT INTEGRATION TESTS ===\n")

required_files <- c(
  "R/pilot_core.R",
  "R/precision_pilot.R",
  "R/lasso_intercept_check.R"
)
if (any(!file.exists(required_files))) {
  stop("Run tests from the project root after unpacking the package.", call. = FALSE)
}

source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/lasso_intercept_check.R")

lasso_intercept_assert_scientific_core()
lasso_intercept_assert_master_infrastructure()
stopifnot(requireNamespace("glmnet", quietly = TRUE))
stopifnot(requireNamespace("batchtools", quietly = TRUE))

beta_test <- c(1, 2, rep(0, 8))
stopifnot(abs(lasso_intercept_signal_variance(beta_test, rho = 0, b = 10L) - 5) < 1e-12)
stopifnot(
  abs(
    lasso_intercept_signal_variance(beta_test, rho = 0.5, b = 10L) -
      (0.5 * 5 + 0.5 * 9)
  ) < 1e-12
)

design_test <- lasso_intercept_prepare_design()
full_settings_test <- lasso_intercept_full_settings(design_test, expected_n_settings = 192L)
selected_test <- lasso_intercept_select_settings(full_settings_test)
tasks_test <- lasso_intercept_make_tasks(
  selected_settings = selected_test,
  R = 50L,
  seed_base = 1500000000L
)

stopifnot(nrow(selected_test) == 12L)
stopifnot(identical(as.integer(selected_test$setting_id), 145:156))
stopifnot(nrow(tasks_test) == 600L)

direct <- lasso_intercept_one_task(
  task_row = tasks_test[1L, , drop = FALSE],
  full_settings = full_settings_test,
  selected_settings = selected_test,
  design = design_test,
  analysis_n_test = 500L
)

stopifnot(nrow(direct) == 1L)
stopifnot(isTRUE(direct$ok[[1L]]))
stopifnot(direct$p[[1L]] == 5000L)
stopifnot(abs(direct$rho[[1L]]) < 1e-12)
stopifnot(direct$n_test[[1L]] == 500L)
stopifnot(direct$n_test_generated[[1L]] >= direct$n_test[[1L]])
stopifnot(is.finite(direct$mse_intercept[[1L]]), direct$mse_intercept[[1L]] > 0)
stopifnot(is.finite(direct$mse_oracle[[1L]]), direct$mse_oracle[[1L]] > 0)
stopifnot(is.finite(direct$mse_lasso[[1L]]), direct$mse_lasso[[1L]] > 0)
stopifnot(is.finite(direct$lambda_min[[1L]]), direct$lambda_min[[1L]] > 0)
stopifnot(is.finite(direct$lambda_1se[[1L]]), direct$lambda_1se[[1L]] > 0)
stopifnot(abs(direct$analytic_snr_rel_error[[1L]]) < 1e-6)

cat("Direct test setting: ", direct$setting_id, " (p=", direct$p, ", rho=", direct$rho, ", SNR=", direct$snr, ")\n", sep = "")
cat("Intercept MSE: ", signif(direct$mse_intercept, 5), "\n", sep = "")
cat("Oracle MSE:    ", signif(direct$mse_oracle, 5), "\n", sep = "")
cat("LASSO MSE:     ", signif(direct$mse_lasso, 5), "\n", sep = "")
cat("LASSO gain:    ", signif(direct$lasso_gain_vs_intercept, 5), "\n", sep = "")
cat("lambda.min:    ", signif(direct$lambda_min, 5), "\n", sep = "")
cat("=== DIRECT INTEGRATION TEST PASSED ===\n\n")

cat("=== BATCHTOOLS / SOCKET SMOKE TEST ===\n")
smoke <- lasso_intercept_batchtools_smoke_test(
  full_settings = full_settings_test,
  selected_settings = selected_test,
  design = design_test,
  task_row = tasks_test[1L, , drop = FALSE],
  analysis_n_test = 500L,
  project_root = ".",
  n_cores = 1L
)
stopifnot(nrow(smoke) == 1L)
stopifnot(isTRUE(smoke$ok[[1L]]))
cat("Batchtools/socket smoke test passed.\n")
cat("=== ALL LASSO / INTERCEPT INTEGRATION TESTS PASSED ===\n\n")

local({
  fault_worker <- lasso_intercept_one_task
  fault_environment <- new.env(parent = environment(fault_worker))
  fault_environment$simulate_data <- function(...) stop("Controlled data-generation failure.")
  environment(fault_worker) <- fault_environment
  fault_task <- tasks_test[1L, , drop = FALSE]
  failed_result <- fault_worker(task_row = fault_task, full_settings = full_settings_test, selected_settings = selected_test, design = design_test, analysis_n_test = 500L)
  stopifnot(isFALSE(failed_result$ok[[1L]]),
            grepl("Controlled data-generation failure", failed_result$error_message[1L], fixed = TRUE))
  for (field in c("task_id", "setting_id", "run_id", "seed")) {
    stopifnot(identical(as.integer(failed_result[[field]][1L]),
                        as.integer(fault_task[[field]][1L])))
  }
  fault_task$seed <- NULL
  malformed_result <- fault_worker(task_row = fault_task, full_settings = full_settings_test, selected_settings = selected_test, design = design_test, analysis_n_test = 500L)
  stopifnot(isFALSE(malformed_result$ok[[1L]]),
            is.na(malformed_result$seed[1L]),
            malformed_result$task_id[1L] == fault_task$task_id[1L],
            malformed_result$setting_id[1L] == fault_task$setting_id[1L],
            malformed_result$run_id[1L] == fault_task$run_id[1L])
})
cat("Structured failure identifiers passed.\n")
