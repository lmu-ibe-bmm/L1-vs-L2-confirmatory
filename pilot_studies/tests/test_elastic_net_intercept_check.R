# Small EN fits with shared data and folds.

cat("\n=== ELASTIC NET / INTERCEPT INTEGRATION TESTS ===\n")

required_files <- c(
  "R/pilot_core.R",
  "R/precision_pilot.R",
  "R/elastic_net_intercept_check.R"
)

if (any(!file.exists(required_files))) {
  stop(
    "Run tests from the project root after unpacking the package.",
    call. = FALSE
  )
}

source("R/pilot_core.R")
source("R/precision_pilot.R")
source("R/elastic_net_intercept_check.R")

elastic_net_intercept_assert_scientific_core()
elastic_net_intercept_assert_master_infrastructure()
stopifnot(requireNamespace("glmnet", quietly = TRUE))
stopifnot(requireNamespace("batchtools", quietly = TRUE))

beta_test <- c(1, 2, rep(0, 8))
stopifnot(
  abs(
    elastic_net_intercept_signal_variance(
      beta_test,
      rho = 0,
      b = 10L
    ) - 5
  ) < 1e-12
)
stopifnot(
  abs(
    elastic_net_intercept_signal_variance(
      beta_test,
      rho = 0.5,
      b = 10L
    ) - (0.5 * 5 + 0.5 * 9)
  ) < 1e-12
)

design_test <- elastic_net_intercept_prepare_design()
full_settings_test <- elastic_net_intercept_full_settings(
  design_test,
  expected_n_settings = 192L
)
selected_test <- elastic_net_intercept_select_settings(
  full_settings_test
)
model_spec_test <- elastic_net_intercept_model_spec()
tasks_test <- elastic_net_intercept_make_tasks(
  selected_settings = selected_test,
  R = 50L,
  seed_base = 1600000000L
)

stopifnot(nrow(selected_test) == 12L)
stopifnot(
  identical(
    as.integer(selected_test$setting_id),
    145:156
  )
)
stopifnot(nrow(tasks_test) == 600L)
stopifnot(nrow(model_spec_test) == 2L)

local({
  fault_worker <- elastic_net_intercept_one_task
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

direct <- elastic_net_intercept_one_task(
  task_row = tasks_test[1L, , drop = FALSE],
  full_settings = full_settings_test,
  selected_settings = selected_test,
  design = design_test,
  model_spec = model_spec_test,
  analysis_n_test = 500L
)

elastic_net_intercept_validate_success_result(
  ans = direct,
  expected_task_id = tasks_test$task_id[1L],
  model_spec = model_spec_test
)

stopifnot(all(direct$p == 5000L))
stopifnot(all(abs(direct$rho) < 1e-12))
stopifnot(all(direct$n_test == 500L))
stopifnot(all(direct$n_test_generated >= direct$n_test))
stopifnot(all(is.finite(direct$mse_intercept)))
stopifnot(all(direct$mse_intercept > 0))
stopifnot(all(is.finite(direct$mse_oracle)))
stopifnot(all(direct$mse_oracle > 0))
stopifnot(all(is.finite(direct$mse_en)))
stopifnot(all(direct$mse_en > 0))
stopifnot(all(is.finite(direct$lambda_min)))
stopifnot(all(direct$lambda_min > 0))
stopifnot(all(is.finite(direct$lambda_1se)))
stopifnot(all(direct$lambda_1se > 0))
stopifnot(all(abs(direct$analytic_snr_rel_error) < 1e-6))

stopifnot(length(unique(direct$data_seed)) == 1L)
stopifnot(length(unique(direct$cv_fold_seed)) == 1L)
stopifnot(max(direct$mse_intercept) - min(direct$mse_intercept) < 1e-15)
stopifnot(max(direct$mse_oracle) - min(direct$mse_oracle) < 1e-15)

cat(
  "Direct test setting: ",
  direct$setting_id[1L],
  " (p=",
  direct$p[1L],
  ", rho=",
  direct$rho[1L],
  ", SNR=",
  direct$snr[1L],
  ")\n",
  sep = ""
)
cat(
  "Intercept MSE: ",
  signif(direct$mse_intercept[1L], 5),
  "\n",
  sep = ""
)
cat(
  "Oracle MSE:    ",
  signif(direct$mse_oracle[1L], 5),
  "\n",
  sep = ""
)

for (i in seq_len(nrow(direct))) {
  cat(
    direct$model[i],
    " alpha=",
    signif(direct$alpha[i], 8),
    ": MSE=",
    signif(direct$mse_en[i], 5),
    ", gain=",
    signif(direct$en_gain_vs_intercept[i], 5),
    ", lambda.min=",
    signif(direct$lambda_min[i], 5),
    "\n",
    sep = ""
  )
}

cat("=== DIRECT INTEGRATION TEST PASSED ===\n\n")

cat("=== BATCHTOOLS / SOCKET SMOKE TEST ===\n")

smoke <- elastic_net_intercept_batchtools_smoke_test(
  full_settings = full_settings_test,
  selected_settings = selected_test,
  design = design_test,
  task_row = tasks_test[1L, , drop = FALSE],
  model_spec = model_spec_test,
  analysis_n_test = 500L,
  project_root = ".",
  n_cores = 1L
)

elastic_net_intercept_validate_success_result(
  ans = smoke,
  expected_task_id = tasks_test$task_id[1L],
  model_spec = model_spec_test
)

cat("Batchtools/socket smoke test passed for both EN alpha values.\n")
cat("=== ALL ELASTIC NET / INTERCEPT INTEGRATION TESTS PASSED ===\n\n")
