###############################################################
# TESTS FOR THE SHARED PILOT CORE
###############################################################
# PURPOSE
# These tests verify the statistical relationships reused by more than one
# pilot: deterministic support construction, population signal-to-noise ratio
# scaling, reproducibility, the ratio-of-means estimand, performance metrics,
# and one small end-to-end model fit.
#
# WHY THESE TESTS MATTER
# An error in this shared layer would propagate into both the precision and
# runtime pilots. The checks therefore focus on relationships that determine
# numerical results, such as support placement, signal-to-noise ratio
# calibration, the paired delta-method Monte Carlo standard error, and the
# exact metrics returned by a fitted model. They do not restate fixed design
# counts or test whether constants are integers.
#
# WHEN TO RUN
# Normally run this file through tests/run_all_tests.R after installing the
# packages with run/00_install_packages.R and before starting any pilot. It can
# also be run alone from the top-level package directory with
#   Rscript tests/test_core.R
# The tests use small examples and never start a full pilot.
###############################################################

source(file.path("R", "pilot_core.R"))

# Print package availability so a skipped model-fitting smoke test is visible.
package_check <- data.frame(
  package = required_packages,
  installed = vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
)
print(package_check)

settings <- make_all_settings(design)

# Use one demanding setting to verify the data-generating mechanism itself
# rather than rechecking every hard-coded design constant.
setting <- settings[
  settings$p == 40 & settings$rho == 0.95 & settings$a == 0.75 &
    settings$k == 5 & settings$snr == 2,
  , drop = FALSE
][1, , drop = FALSE]

beta_1 <- make_beta(setting, seed = 123456L, design = design)
beta_2 <- make_beta(setting, seed = 123456L, design = design)
beta_3 <- make_beta(setting, seed = 123457L, design = design)

# Equal seeds reproduce coefficients exactly. A different seed changes
# magnitudes or signs but not the deterministic support.
stopifnot(
  identical(beta_1$beta0, beta_2$beta0),
  identical(beta_1$beta, beta_2$beta),
  identical(beta_1$block_signs, beta_2$block_signs),
  identical(beta_1$active_blocks, beta_3$active_blocks),
  identical(beta_1$active_variables, beta_3$active_variables),
  !identical(beta_1$beta0, beta_3$beta0),
  identical(beta_1$active_blocks, seq_len(setting$B_active)),
  abs(beta_1$achieved_snr_theoretical - setting$snr) < 1e-10
)

# Within every active block, the first k variables are active and share one
# sign. This is the central structural feature of the coefficient-generating
# mechanism.
for (block_index in seq_along(beta_1$active_blocks)) {
  block <- beta_1$active_blocks[block_index]
  block_columns <- ((block - 1L) * setting$b + 1L):(block * setting$b)
  expected_active <- block_columns[seq_len(setting$k)]
  observed_active <- intersect(block_columns, beta_1$active_variables)
  raw_coefficients <- beta_1$beta0[expected_active]

  stopifnot(
    identical(expected_active, observed_active),
    all(raw_coefficients > 0) || all(raw_coefficients < 0),
    all(abs(raw_coefficients) >= design$beta_magnitude_min),
    all(abs(raw_coefficients) <= design$beta_magnitude_max)
  )
}

# The primary effect is a ratio of means, not a mean of repetition-specific
# ratios. The unequal prediction mean squared error pairs make the difference
# observable.
lasso_mse <- c(1, 10)
ridge_mse <- c(0.9, 9.5)
ratio_summary <- ratio_of_means_summary(lasso_mse, ridge_mse)
mean_of_ratios <- mean((lasso_mse - ridge_mse) / lasso_mse)

stopifnot(
  abs(
    ratio_summary$theta -
      ((mean(lasso_mse) - mean(ridge_mse)) / mean(lasso_mse))
  ) < 1e-14,
  abs(ratio_summary$theta - mean_of_ratios) > 1e-3,
  is.finite(ratio_summary$tau_influence),
  is.finite(ratio_summary$mcse)
)

# Compare the influence-value Monte Carlo standard error with the direct paired
# delta-method form.
mean_lasso <- mean(lasso_mse)
mean_ridge <- mean(ridge_mse)
delta_variance <- (
  (mean_ridge^2 / mean_lasso^4) * var(lasso_mse) +
    (1 / mean_lasso^2) * var(ridge_mse) -
    2 * (mean_ridge / mean_lasso^3) * cov(lasso_mse, ridge_mse)
) / length(lasso_mse)
stopifnot(abs(ratio_summary$mcse^2 - delta_variance) < 1e-12)

# Check the small metric helpers with values whose answers are available by
# hand. These functions are used in every fitted-model task.
prediction_metrics <- compute_prediction_metrics(
  y = c(1, 2, 3),
  prediction = c(1, 1, 4)
)
stopifnot(
  abs(prediction_metrics$mse - 2 / 3) < 1e-14,
  abs(prediction_metrics$mae - 2 / 3) < 1e-14
)

coefficient_metrics <- compute_coefficient_mse(
  beta_hat = c(1, 0.5, 0),
  beta_true = c(1, 1, 0),
  active_variables = c(1L, 2L)
)
stopifnot(
  abs(coefficient_metrics$coefficient_mse_overall - 0.25 / 3) < 1e-14,
  abs(coefficient_metrics$coefficient_mse_active - 0.125) < 1e-14,
  abs(coefficient_metrics$coefficient_mse_inactive) < 1e-14
)

# If all model packages are installed, perform one small real fit with shared
# folds. Reducing n_test keeps this verification fast without changing the
# data-generating mechanism.
if (all(package_check$installed)) {
  smoke_setting <- settings[
    settings$p == 40 & settings$rho == 0.6,
    , drop = FALSE
  ][1, , drop = FALSE]
  smoke_setting$n_test <- 100L

  beta_object <- make_beta(smoke_setting, 222001L, design)
  data_object <- simulate_data(smoke_setting, beta_object, 222002L, design)
  fold_id <- make_cv_foldid(
    nrow(data_object$X_train), design$glmnet_nfolds, 222003L
  )

  fitted_models <- list(
    ridge = fit_glmnet(
      data_object, 0, fold_id,
      beta_object$beta, beta_object$active_variables, design
    ),
    lasso = fit_glmnet(
      data_object, 1, fold_id,
      beta_object$beta, beta_object$active_variables, design
    ),
    elastic_net = fit_glmnet(
      data_object, design$en_alpha, fold_id,
      beta_object$beta, beta_object$active_variables, design
    ),
    elastic_net_095 = fit_glmnet(
      data_object, design$en_near_lasso_alpha, fold_id,
      beta_object$beta, beta_object$active_variables, design
    )
  )

  for (fitted_model in fitted_models) {
    stopifnot(
      is.finite(fitted_model$mse),
      is.finite(fitted_model$mae),
      is.finite(fitted_model$coefficient_mse_overall),
      is.finite(fitted_model$coefficient_mse_active),
      is.finite(fitted_model$coefficient_mse_inactive),
      length(fitted_model$beta_hat) == smoke_setting$p
    )
  }

  random_forest <- fit_rf(data_object, smoke_setting$p, 222004L, design)
  stopifnot(is.finite(random_forest$mse), is.finite(random_forest$mae))
}

cat("All shared-core tests passed.\n")
