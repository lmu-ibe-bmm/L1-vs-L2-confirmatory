###############################################################
# TESTS FOR THE BLINDED PRECISION PILOT
###############################################################
# PURPOSE
# These tests verify the variance-only Random Forest surrogate, the derivation
# of the influence-value standard deviation and the required number of
# repetitions, the exact data transformation used in the report figure, both
# figure-source modes, one real precision task, and the small serial
# execution path.
#
# WHY THESE TESTS MATTER
# The pilot is intentionally blinded, so its validity depends on keeping the
# raw Random-Forest mean out of the planning effect while preserving the paired
# dependence pattern. The tests use hand-checkable vectors to verify centering,
# rescaling, the assumed mean shift, and the conversion from influence-value
# standard deviation to a repetition requirement. Synthetic figure data
# verify ratio direction and threshold counts. Each selected input receives
# its own output name.
#
# WHEN TO RUN
# Run through tests/run_all_tests.R after package installation and before
# run/02_run_precision_pilot.R. To run only this file from the top-level package
# directory:
#   Rscript tests/test_precision_pilot.R
# No full precision pilot is started here.
###############################################################

source(file.path("R", "precision_pilot.R"))

# Test timestamp-based selection and explicit input paths.
selection_directory <- tempfile("precision_source_selection_")
dir.create(selection_directory)
first_summary <- file.path(
  selection_directory,
  "sd_precision_setting_summary_2000-01-01_00-00-00.csv"
)
second_summary <- file.path(
  selection_directory,
  "sd_precision_setting_summary_2000-01-02_00-00-00.csv"
)
file.create(first_summary, second_summary)

latest_source <- resolve_precision_figure_source(
  "latest",
  results_dir = selection_directory
)
explicit_source <- resolve_precision_figure_source(first_summary)

stopifnot(
  identical(latest_source$mode, "latest"),
  identical(basename(latest_source$input_file), basename(second_summary)),
  identical(
    latest_source$figure_file_stem,
    "Supp_Fig_RF_LASSO_SD_ratio_boxplot_pretty_2000-01-02_00-00-00"
  ),
  identical(explicit_source$mode, "file"),
  identical(basename(explicit_source$input_file), basename(first_summary)),
  identical(
    explicit_source$figure_file_stem,
    "Supp_Fig_RF_LASSO_SD_ratio_boxplot_pretty_2000-01-01_00-00-00"
  )
)
unlink(selection_directory, recursive = TRUE)

# The surrogate mean is controlled only by theta_assumed. Centered Random Forest
# deviations determine the correlation pattern but not the mean level.
lasso_mse <- c(1.00, 1.10, 0.90, 1.05, 0.95, 1.02)
rf_mse <- c(2.00, 2.30, 1.70, 2.10, 1.90, 2.20)

surrogate_zero <- make_variance_only_rf_surrogate(
  lasso_mse, rf_mse, sd_multiplier = 1.0, theta_assumed = 0
)
surrogate_positive <- make_variance_only_rf_surrogate(
  lasso_mse, rf_mse, sd_multiplier = 1.0, theta_assumed = 0.02
)
surrogate_negative <- make_variance_only_rf_surrogate(
  lasso_mse, rf_mse, sd_multiplier = 1.0, theta_assumed = -0.02
)

stopifnot(
  abs(mean(surrogate_zero$mse_surrogate) - mean(lasso_mse)) < 1e-12,
  abs(mean(surrogate_positive$mse_surrogate) - 0.98 * mean(lasso_mse)) < 1e-12,
  abs(mean(surrogate_negative$mse_surrogate) - 1.02 * mean(lasso_mse)) < 1e-12,
  abs(sd(surrogate_zero$mse_surrogate) - sd(lasso_mse)) < 1e-12,
  abs(
    cor(lasso_mse, surrogate_zero$mse_surrogate) -
      cor(lasso_mse, rf_mse)
  ) < 1e-12
)

# Verify that the Random Forest deviations can also be scaled upward when the
# raw Random Forest pilot variability is smaller than the LASSO variability.
low_variability_rf <- 1 + 0.2 * (lasso_mse - mean(lasso_mse))
surrogates <- lapply(c(1.0, 1.1, 1.2, 1.3), function(multiplier) {
  return(make_variance_only_rf_surrogate(
    lasso_mse,
    low_variability_rf,
    sd_multiplier = multiplier,
    theta_assumed = 0
  ))
})

for (index in seq_along(surrogates)) {
  multiplier <- c(1.0, 1.1, 1.2, 1.3)[index]
  stopifnot(
    abs(
      sd(surrogates[[index]]$mse_surrogate) -
        multiplier * sd(lasso_mse)
    ) < 1e-12
  )
}
stopifnot(all(diff(vapply(surrogates, `[[`, numeric(1), "scale_factor")) > 0))

# Synthetic results exercise the complete setting-level planning summary
# without fitting models.
small_settings <- make_all_settings(design)[1:2, , drop = FALSE]
set.seed(42)
synthetic_results <- bind_rows(lapply(1:2, function(setting_id) {
  repetitions <- 100L
  lasso <- 1 + 0.1 * rnorm(repetitions)
  rf <- 1.5 + 0.14 * scale(lasso)[, 1] + 0.05 * rnorm(repetitions)

  return(tibble(
    task_id = as.integer((setting_id - 1L) * repetitions + seq_len(repetitions)),
    setting_id = setting_id,
    run_id = seq_len(repetitions),
    ok = TRUE,
    mse_lasso = lasso,
    mse_rf = rf
  ))
}))

precision_summary <- summarise_precision_pilot(
  synthetic_results,
  small_settings,
  design,
  sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3),
  theta_assumed_grid = c(-0.02, 0, 0.02)
)

stopifnot(
  all(is.finite(precision_summary$setting_summary$surrogate_tau_max)),
  all(
    precision_summary$setting_summary$R_required_sensitivity >= design$R_min
  ),
  is.finite(precision_summary$global_summary$R_global_sensitivity),
  identical(
    precision_summary$sensitivity_summary$sd_multiplier,
    c(1.0, 1.1, 1.2, 1.3)
  ),
  all(diff(
    precision_summary$sensitivity_summary$R_global_unrounded
  ) >= 0),
  max(precision_summary$sensitivity_summary$R_global_unrounded) ==
    precision_summary$global_summary$R_global_sensitivity,
  identical(
    precision_summary$global_summary$sd_multiplier_grid,
    "1,1.1,1.2,1.3"
  )
)

# Hand-checkable data verify ratio direction and threshold counts.
figure_input <- data.frame(
  setting_id = 1:4, sd_mse_lasso = c(1, 2, 1, 1),
  sd_mse_rf_raw = c(0.5, 2, 1.5, 3)
)
figure_prepared <- prepare_rf_lasso_sd_ratio(figure_input)
stopifnot(
  abs(figure_prepared$summary$median_ratio - 1.25) < 1e-12,
  figure_prepared$summary$n_above_1 == 2L,
  figure_prepared$summary$n_above_1_3 == 2L,
  figure_prepared$summary$n_above_2 == 1L
)

# Check Portable Network Graphics and Portable Document Format export when the
# plotting package is available. The numeric plot transformation above is
# tested independently of package availability.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  figure_directory <- tempfile("figure_test_")
  figure_files <- save_rf_lasso_sd_ratio_figure(
    figure_input,
    output_dir = figure_directory,
    file_stem = "rf_lasso_sd_test"
  )
  stopifnot(
    file.exists(figure_files$png),
    file.exists(figure_files$pdf),
    file.info(figure_files$png)$size > 0,
    file.info(figure_files$pdf)$size > 0
  )
}

packages_available <- all(vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
))

# One real task verifies the connection among the data-generating mechanism,
# LASSO, and Random Forest.
if (packages_available) {
  smoke_setting <- make_all_settings(design)
  smoke_setting <- smoke_setting[
    smoke_setting$p == 40 & smoke_setting$rho == 0.95 &
      smoke_setting$a == 0.75 & smoke_setting$k == 10,
    , drop = FALSE
  ][1, , drop = FALSE]
  smoke_setting$n_test <- 100L
  smoke_setting$task_id <- 1L
  smoke_setting$run_id <- 1L
  smoke_setting$seed <- 880001L

  task_result <- run_one_sd_task(smoke_setting, design)
  stopifnot(
    isTRUE(task_result$ok),
    is.finite(task_result$mse_lasso),
    is.finite(task_result$mse_rf),
    !any(c("mae_lasso", "mae_rf") %in% names(task_result))
  )

  # Exercise the user-facing serial orchestration path with one setting and two
  # repetitions. The task-level checks above already validate the values.
  tiny_design <- design
  tiny_design$n_values <- 200
  tiny_design$p_over_n_values <- 0.2
  tiny_design$rho_values <- 0
  tiny_design$a_values <- 0.25
  tiny_design$k_values <- 1
  tiny_design$snr_values <- 0.5

  small_pilot <- run_precision_pilot(
    R_sd = 2L,
    backend = "serial",
    design = tiny_design,
    save_outputs = FALSE,
    sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3)
  )
  stopifnot(
    small_pilot$global_summary$n_success == 2L,
    !"run_results" %in% names(small_pilot),
    !"tasks" %in% names(small_pilot)
  )
}

cat("All precision-pilot tests passed.\n")
