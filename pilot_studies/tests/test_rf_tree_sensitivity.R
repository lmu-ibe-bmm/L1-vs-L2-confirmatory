# Small paired-tree and planning checks.

source(file.path("R", "rf_tree_sensitivity.R"))

stopifnot(
  is.function(run_one_rf_tree_sensitivity_task),
  is.function(summarise_rf_tree_sensitivity),
  is.function(save_rf_tree_sensitivity_figure),
  is.function(run_rf_tree_sensitivity)
)

test_design <- make_design()
test_design$rf_tree_counts <- c(200L, 500L, 1000L)

test_settings <- make_all_settings(test_design)[1:2, , drop = FALSE]
repetitions <- 8L

synthetic_rows <- vector("list", nrow(test_settings) * repetitions)
row_index <- 1L

for (setting_index in seq_len(nrow(test_settings))) {
  setting_id <- test_settings$setting_id[setting_index]

  for (run_id in seq_len(repetitions)) {

    alternating <- if (run_id %% 2L == 0L) 1 else -1
    lasso_pmse <- 1 + 0.015 * run_id + 0.01 * setting_index

    synthetic_rows[[row_index]] <- tibble::tibble(
      task_id = row_index,
      setting_id = setting_id,
      run_id = run_id,
      ok = TRUE,
      mse_lasso = lasso_pmse,
      mse_rf_200 = 1.02 + 1.18 * (lasso_pmse - 1) + 0.020 * alternating,
      mse_rf_500 = 1.02 + 1.10 * (lasso_pmse - 1) + 0.012 * alternating,
      mse_rf_1000 = 1.02 + 1.06 * (lasso_pmse - 1) + 0.008 * alternating
    )
    row_index <- row_index + 1L
  }
}

synthetic_results <- dplyr::bind_rows(synthetic_rows)
synthetic_summary <- summarise_rf_tree_sensitivity(
  run_results = synthetic_results,
  settings = test_settings,
  design = test_design
)

stopifnot(
  nrow(synthetic_summary$setting_summary) == 6L,
  nrow(synthetic_summary$tree_count_summary) == 3L,
  nrow(synthetic_summary$sensitivity_summary) == 12L,
  nrow(synthetic_summary$global_summary) == 3L,
  nrow(synthetic_summary$comparison_with_200) == 6L,
  identical(
    sort(unique(synthetic_summary$setting_summary$num_trees)),
    c(200L, 500L, 1000L)
  ),
  identical(
    synthetic_summary$tree_count_summary$num_trees,
    c(200L, 500L, 1000L)
  ),
  all(is.finite(
    synthetic_summary$tree_count_summary$median_sd_ratio_rf_lasso
  )),
  all(is.finite(
    synthetic_summary$tree_count_summary$median_correlation_lasso_rf
  )),
  all(c(
    "change_sd_ratio_from_200",
    "change_correlation_from_200",
    "change_sensitivity_R_from_200"
  ) %in% names(synthetic_summary$comparison_with_200))
)

baseline_rows <- synthetic_summary$comparison_with_200$num_trees == 200L
stopifnot(
  all(
    synthetic_summary$comparison_with_200$change_sd_mse_rf_from_200[
      baseline_rows
    ] == 0
  ),
  all(
    synthetic_summary$comparison_with_200$change_sd_ratio_from_200[
      baseline_rows
    ] == 0
  ),
  all(
    synthetic_summary$comparison_with_200$change_correlation_from_200[
      baseline_rows
    ] == 0
  ),
  all(
    synthetic_summary$comparison_with_200$change_raw_R_from_200[
      baseline_rows
    ] == 0
  ),
  all(
    synthetic_summary$comparison_with_200$change_sensitivity_R_from_200[
      baseline_rows
    ] == 0
  )
)

public_names <- unique(c(
  names(synthetic_summary$setting_summary),
  names(synthetic_summary$tree_count_summary),
  names(synthetic_summary$sensitivity_summary),
  names(synthetic_summary$global_summary),
  names(synthetic_summary$comparison_with_200)
))
stopifnot(
  !any(c(
    "mse_lasso", "mse_rf_200", "mse_rf_500", "mse_rf_1000"
  ) %in% public_names),
  all(c(
    "sd_mse_lasso", "sd_mse_rf_raw", "corr_lasso_rf_raw"
  ) %in% names(synthetic_summary$setting_summary))
)

temporary_figure_directory <- tempfile("rf_tree_sensitivity_test_")
dir.create(temporary_figure_directory)

test_figure <- save_rf_tree_sensitivity_figure(
  setting_summary = synthetic_summary$setting_summary,
  output_dir = temporary_figure_directory,
  file_stem = "rf_tree_sensitivity_test"
)

stopifnot(
  file.exists(test_figure$png),
  file.exists(test_figure$pdf),
  file.info(test_figure$png)$size > 0,
  file.info(test_figure$pdf)$size > 0
)
unlink(temporary_figure_directory, recursive = TRUE)

smoke_design <- make_design()
smoke_design$n_test <- 40L
smoke_design$rf_tree_counts <- c(10L, 20L)

smoke_setting <- make_all_settings(smoke_design)
smoke_setting <- smoke_setting[
  smoke_setting$p == 40 &
    smoke_setting$rho == 0.6 &
    smoke_setting$a == 0.25 &
    smoke_setting$k == 1 &
    smoke_setting$snr == 0.5,
  ,
  drop = FALSE
][1, , drop = FALSE]

smoke_task <- make_tasks(
  settings = smoke_setting,
  R = 1L,
  seed_base = 9200L
)
smoke_result <- run_one_rf_tree_sensitivity_task(
  task_row = smoke_task,
  design = smoke_design
)

stopifnot(
  isTRUE(smoke_result$ok),
  is.na(smoke_result$failed_step),
  is.finite(smoke_result$mse_lasso),
  is.finite(smoke_result$mse_rf_10),
  is.finite(smoke_result$mse_rf_20),
  smoke_result$mtry == floor(sqrt(smoke_setting$p)),
  smoke_result$time_rf_fit_10_sec >= 0,
  smoke_result$time_rf_fit_20_sec >= 0
)

cat("Random Forest tree-count sensitivity tests passed.\n")
