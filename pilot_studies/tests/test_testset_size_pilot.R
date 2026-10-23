###############################################################
# TESTS FOR THE TEST-SET-SIZE PILOT
###############################################################
# PURPOSE
# These tests ensure that the pilot's compact local data-generating mechanism
# uses the same structural design as the shared core. They also verify
# reproducibility and support construction, validate the exact Gaussian
# test-error formulas, and exercise one small glmnet evaluation.
#
# WHY THESE TESTS MATTER
# The test-set pilot intentionally has a small standalone data-generating
# mechanism, so the highest risk is silent drift from the shared covariance,
# support, or signal-to-noise ratio definitions.
# The first checks compare those relationships directly. The remaining checks
# verify the analytical variance formulas with a perfect fitted model. A small
# real fit then verifies the nested prediction path without rerunning either
# complete pilot component.
#
# WHEN TO RUN
# Run through tests/run_all_tests.R after package installation and before
# run/01_run_testset_size_pilot.R. To run only this file from the top-level
# package directory:
#   Rscript tests/test_testset_size_pilot.R
# The examples are deliberately small and do not start the complete pilot.
###############################################################

# Load the shared core first so the two data-generating implementations can be
# compared.
source(file.path("R", "pilot_core.R"))
source(file.path("R", "testset_size_pilot.R"))

testset_design <- make_testset_pilot_design()
settings <- make_testset_settings(testset_design)

setting <- settings[
  settings$p == 40 & settings$rho == 0.95 &
    settings$a == 0.75 & settings$k == 5 & settings$snr == 2,
  , drop = FALSE
][1, , drop = FALSE]

testset_beta_1 <- make_beta_testset(
  setting,
  seed = 12345,
  design = testset_design
)
testset_beta_2 <- make_beta_testset(
  setting,
  seed = 12345,
  design = testset_design
)

# Compare the derived setting and deterministic support with the shared core.
# The two implementations deliberately manage random draws independently, so
# equality of seed-specific magnitudes or signs is not a scientific requirement.
core_design <- make_design()
core_setting <- make_setting(
  n = setting$n,
  p_over_n = setting$p_over_n,
  rho = setting$rho,
  a = setting$a,
  k = setting$k,
  snr = setting$snr,
  design = core_design
)
core_beta <- make_beta(core_setting, seed = 12345, design = core_design)
stopifnot(
  core_setting$p == setting$p,
  core_setting$b == setting$b,
  core_setting$G == setting$G,
  core_setting$B_active == setting$B_active,
  core_setting$q_active == setting$q_active,
  identical(core_beta$active_blocks, testset_beta_1$active_blocks),
  identical(core_beta$active_variables, testset_beta_1$active_variables)
)

# Equal seeds reproduce all random coefficient components exactly.
stopifnot(
  identical(testset_beta_1$active_blocks, testset_beta_2$active_blocks),
  identical(testset_beta_1$active_variables, testset_beta_2$active_variables),
  identical(testset_beta_1$beta0, testset_beta_2$beta0),
  identical(testset_beta_1$block_signs, testset_beta_2$block_signs),
  identical(testset_beta_1$active_blocks, seq_len(setting$B_active))
)

# The first k positions are active in every active block and share one sign.
for (block_index in seq_along(testset_beta_1$active_blocks)) {
  block <- testset_beta_1$active_blocks[block_index]
  block_columns <- ((block - 1L) * setting$b + 1L):(block * setting$b)
  expected_active <- block_columns[seq_len(setting$k)]
  observed_active <- intersect(
    block_columns,
    testset_beta_1$active_variables
  )
  raw_coefficients <- testset_beta_1$beta0[expected_active]

  stopifnot(
    identical(observed_active, expected_active),
    all(raw_coefficients > 0) || all(raw_coefficients < 0),
    all(abs(raw_coefficients) >= testset_design$beta_magnitude_min),
    all(abs(raw_coefficients) <= testset_design$beta_magnitude_max)
  )
}

achieved_snr <- quadform_block_equicorr_testset(
  testset_beta_1$beta,
  setting$b,
  setting$rho
) / testset_design$sigma2
stopifnot(abs(achieved_snr - setting$snr) < 1e-10)

# Verify the theoretical inverse-n variance and inverse-square-root-n
# standard-deviation laws.
scaling <- make_theoretical_n_test_scaling_table()
stopifnot(
  abs(scaling$var_ratio_vs_5000[scaling$n_test == 5000] - 1) < 1e-15,
  abs(scaling$var_ratio_vs_5000[scaling$n_test == 10000] - 0.5) < 1e-15,
  abs(scaling$var_ratio_vs_5000[scaling$n_test == 20000] - 0.25) < 1e-15,
  abs(
    scaling$sd_ratio_vs_5000[scaling$n_test == 10000] - sqrt(0.5)
  ) < 1e-15
)

# With a perfect fitted regression function, the prediction error is exactly
# normally distributed with mean zero and variance one, so its expected mean
# squared error is one and the variance of that mean squared error is
# 2 / n_test.
perfect_fit <- list(
  beta_hat = testset_beta_1$beta,
  intercept_hat = 0
)
moments <- conditional_test_error_moments_testset(
  beta_true = testset_beta_1$beta,
  fit_obj = perfect_fit,
  b = setting$b,
  rho = setting$rho,
  n_test_values = c(5000L, 10000L, 20000L),
  sigma2 = 1
)
stopifnot(
  all(abs(moments$expected_mse - 1) < 1e-12),
  all(abs(moments$var_test_mse - 2 / moments$n_test) < 1e-12),
  all(abs(moments$cv_test_mse - sqrt(2 / moments$n_test)) < 1e-12)
)

# If glmnet is installed, fit all four models once and evaluate two small
# nested test-set sizes. This checks the real prediction path without running
# the full exact or runtime arm.
if (requireNamespace("glmnet", quietly = TRUE)) {
  smoke_setting <- settings[
    settings$p == 40 & settings$rho == 0.6,
    , drop = FALSE
  ][1, , drop = FALSE]
  fixed_models <- fit_fixed_models_testset(
    setting_row = smoke_setting,
    base_seed = 777000L,
    design = testset_design
  )

  for (method in names(fixed_models$fits)) {
    stopifnot(
      length(fixed_models$fits[[method]]$beta_hat) == smoke_setting$p,
      is.finite(fixed_models$fits[[method]]$intercept_hat),
      is.finite(fixed_models$fits[[method]]$lambda_min)
    )
  }

  small_design <- testset_design
  small_design$n_test_values <- c(100L, 200L)
  small_design$runtime_chunk_size <- 100L
  evaluation <- run_nested_test_evaluation_testset(
    setting_row = smoke_setting,
    fixed = fixed_models,
    test_seed = 778000L,
    design = small_design
  )
  stopifnot(
    all(is.finite(evaluation$mse)),
    all(is.finite(evaluation$mae)),
    all(evaluation$mse >= 0),
    all(evaluation$mae >= 0)
  )
}

cat("All test-set-size-pilot tests passed.\n")
