###############################################################
# STANDALONE TEST-SET-SIZE PILOT
#
# PURPOSE
# This script quantifies only the additional variation in prediction mean
# squared error and prediction mean absolute error caused by a finite test set
# while the fitted model is held fixed. It covers Ridge (alpha = 0), Elastic
# Net (1/3), near-LASSO Elastic Net
# (0.95 = 1 - epsilon, epsilon = 0.05), and LASSO (1).
# In later comments and column names, MSE means mean squared error, MAE means
# mean absolute error, SD means standard deviation, and SNR means
# signal-to-noise ratio.
#
# CONTRIBUTION TO THE OVERALL STUDY
# This pilot answers one narrow design question, whether a common independent
# test-set size of 5,000 makes the reported prediction errors sufficiently
# stable. It does not choose the number of Monte Carlo repetitions and it does
# not use method differences to investigate confirmatory claims. Holding each
# fitted model fixed isolates test-sample noise from training-sample noise and
# from the between-repetition uncertainty studied by the precision pilot.
#
# WHY A SMALL LOCAL DATA-GENERATING MECHANISM IS USED
# The pilot holds a trained model fixed and then studies only the additional
# test-set variation. It therefore remains runnable without batchtools and
# ranger. Its compact, separately named data-generating helpers are kept in
# this file. tests/test_testset_size_pilot.R compares their output with the
# shared core so that the specialized implementation stays aligned.
#
# TWO COMPLEMENTARY ARMS
#   A) An exact variance arm. After one fixed model fit, the expectation and
#      variance of prediction mean squared error and prediction mean absolute
#      error are derived analytically from the Gaussian
#      prediction-error distribution.
#   B) A matrix/prediction runtime arm for p = 40 and p = 5,000. The same fit
#      is evaluated on newly simulated test data generated in memory-efficient
#      chunks.
#
# STATISTICAL BASIS AND LITERATURE
# Conditional on a fitted model, prediction mean squared error is the average
# of independent squared prediction errors. Its variance is therefore
# Var(e^2) / n_test. For Gaussian e with mean m and variance v,
# E(e^2) = v + m^2 and Var(e^2) = 2 v^2 + 4 m^2 v. The prediction mean absolute
# error formulas below use the corresponding folded-normal moments. This
# analytical arm avoids adding another simulation layer when the conditional
# moments are available exactly.
#
# The separation of aims, data-generating mechanisms, estimands, methods, and
# performance measures follows Morris et al. (2019) and the planning and
# reporting framework of Siepe et al. (2024). The conditional variance formulas
# above follow directly from averaging independent test losses and from
# Gaussian moments. These references do not prescribe the four candidate
# test-set sizes, the use of one fixed fit in the exact component, the five
# runtime evaluations, or the final choice of 5,000 observations. Those are
# study-specific design decisions.
#
# HOW TO INTERPRET THE OUTPUT
# The exact component describes the conditional standard deviation and
# coefficient of variation of prediction mean squared error and prediction
# mean absolute error for each candidate size. The runtime component shows the
# computing cost of generating test data and predictions at those sizes. The
# public files do not contain individual prediction-error values. The pilot
# does not perform a hypothesis test or apply an automatic pass/fail threshold.
# The choice of n_test = 5,000 is a documented precision-cost decision based on
# the two complementary summaries.
#
# WHEN AND HOW TO RUN THIS SCRIPT
# From the top-level package directory, first run run/00_install_packages.R, then
# tests/run_all_tests.R, and finally run/01_run_testset_size_pilot.R. The
# launcher writes variability tables, runtime tables, theoretical scaling, and
# one summary R data object to outputs/testset_size_pilot/. All files from one
# execution share one timestamp, which also identifies its software-session
# record. Sourcing this file defines the functions but does not start the pilot
# automatically.
###############################################################

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

required_packages_testset <- "glmnet"

# Locate the base-R helper through the actual module path. The standalone pilot
# still needs only glmnet, and can be sourced by an absolute path from elsewhere.
testset_source_files <- vapply(sys.frames(), function(frame) {
  return(if (!is.null(frame$ofile)) as.character(frame$ofile) else NA_character_)
}, character(1))
testset_module_files <- testset_source_files[
  !is.na(testset_source_files) &
    basename(testset_source_files) == "testset_size_pilot.R"
]
if (length(testset_module_files) == 0L) {
  testset_command_file <- sub("--file=", "", commandArgs()[
    grepl("^--file=", commandArgs())
  ])
  testset_module_files <- c(testset_command_file, "R/testset_size_pilot.R")
  testset_module_files <- testset_module_files[file.exists(testset_module_files)]
}
if (length(testset_module_files) == 0L) {
  stop("Cannot locate testset_size_pilot.R and its reproducibility helper.")
}
source(file.path(dirname(tail(testset_module_files, 1L)), "reproducibility.R"))

# PURPOSE
# Verify that the only external package required by this standalone pilot is
# available before a long calculation starts.
#
# WHY THIS FUNCTION EXISTS
# A missing glmnet installation should produce one direct message at the start
# of a run, not an obscure failure after settings and data have been created.
# Different installed versions are reported as warnings and do not stop a run.
#
# INPUTS
# None. The package names are stored in required_packages_testset above.
#
# RETURN VALUE
# invisible(TRUE) when every package is available, otherwise the function
# stops and lists the missing package names. Reference-version differences warn.
check_required_packages_testset <- function() {
  # The shared base-R checker keeps missing-package errors separate from
  # differences in installed R, glmnet, and its Matrix dependency.
  check_reference_environment(packages = c(required_packages_testset, "Matrix"))
  return(invisible(TRUE))
}

# PURPOSE
# Create the complete set of fixed constants used by the test-set-size pilot.
#
# WHY THIS FUNCTION EXISTS
# Keeping factor levels, data-generating parameters, model settings, candidate
# test-set sizes, and seeds in one named object makes the pilot reproducible and
# keeps assumptions visible at every function boundary.
#
# INPUTS
# None. This package documents one fixed pilot design.
#
# RETURN VALUE
# A named list used by all functions below. It contains the two dimensional
# extremes, four candidate n_test values, glmnet settings, and separate seeds
# and repetition counts for the exact and runtime arms.
make_testset_pilot_design <- function() {
  return(list(
    design_name = "ridge_lasso_gaussian_linear_testset_size_pilot",
    n_values = c(200L, 1000L),
    p_over_n_values = c(0.2, 5),
    rho_values = c(0, 0.6, 0.8, 0.95),
    a_values = c(0.25, 0.75),
    k_values = c(1L, 5L, 10L),
    snr_values = c(0.5, 2),
    b = 10L,
    sigma2 = 1,
    beta_magnitude_min = 0.5,
    beta_magnitude_max = 1.5,
    glmnet_nfolds = 5L,
    glmnet_type_measure = "mse",
    glmnet_standardize = TRUE,
    glmnet_intercept = TRUE,
    glmnet_prediction_s = "lambda.min",
    # Equal coefficients on the L1 and squared-L2 penalty terms in glmnet.
    en_alpha = 1 / 3,
    en_near_lasso_alpha = 0.95,
    # The comparison focuses on the two dimensional extremes. This fixed
    # project choice covers the smallest and largest predictor matrices while
    # retaining every structural factor level within each extreme.
    p_values_for_testset_pilot = c(40L, 5000L),
    # The four candidate sizes are prespecified project choices. The smallest
    # value provides a lower reference, while 5,000, 10,000, and 20,000 show
    # the precision and runtime changes from successively larger evaluations.
    n_test_values = c(2500L, 5000L, 10000L, 20000L),
    runtime_chunk_size = 2500L,
    # Exact arm: one fixed training fit per setting is enough for the intended
    # conditional question, because the fitted model is deliberately held
    # fixed while finite-test-set uncertainty is calculated analytically.
    R_base_exact = 1L,
    # Runtime benchmark: five evaluations describe test-data generation and
    # prediction cost. This arm does not estimate between-training-sample
    # variation and therefore does not require a study-sized repetition count.
    R_runtime_eval = 5L,
    seed_base_exact = 500000000L,
    seed_base_runtime = 600000000L,
    seed_base_evaluation = 650000000L
  ))
}

# PURPOSE
# Expand the test-set pilot's factor grid and retain only the two specified
# dimensional extremes, p = 40 and p = 5,000.
#
# WHY THIS FUNCTION EXISTS
# The test-set question does not require every n-by-p/n combination from the
# complete factorial design. Studying the smallest and largest predictor matrices
# covers the two relevant computational extremes while all correlation,
# sparsity, and signal levels remain represented.
#
# INPUT
# design: The list returned by make_testset_pilot_design().
#
# RETURN VALUE
# A data frame with one row per retained setting. Besides the experimental
# factors it includes p, number of blocks, active blocks, active variables,
# realized sparsity, and a stable setting_id. The derived columns are retained
# in the saved pilot object so each result can be interpreted without manually
# reconstructing its support size.
make_testset_settings <- function(design = make_testset_pilot_design()) {
  # Step 1: combine n with p/n and calculate the implied predictor dimension.
  np <- expand.grid(
    n = design$n_values,
    p_over_n = design$p_over_n_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  np$p <- as.integer(round(np$n * np$p_over_n))
  np <- np[np$p %in% design$p_values_for_testset_pilot, , drop = FALSE]

  # Step 2: cross the retained dimensions with all structural design factors.
  settings <- merge(
    np,
    expand.grid(
      rho = design$rho_values,
      a = design$a_values,
      k = design$k_values,
      snr = design$snr_values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    by = NULL
  )

  # Step 3: derive support sizes from block size b, active-block share a, and
  # active variables per active block k.
  settings$b <- as.integer(design$b)
  settings$G <- as.integer(settings$p / settings$b)
  settings$B_active <- as.integer(round(settings$a * settings$G))
  settings$q_active <- as.integer(settings$B_active * settings$k)
  settings$s_realized <- settings$q_active / settings$p

  # Step 4: impose a deterministic order before assigning identifiers. This
  # keeps seeds and output rows stable across machines.
  settings <- settings[order(
    settings$p,
    settings$rho,
    settings$a,
    settings$k,
    settings$snr
  ), , drop = FALSE]
  rownames(settings) <- NULL
  settings$setting_id <- seq_len(nrow(settings))

  return(settings)
}

# PURPOSE
# Assign every training observation to one of the shared cross-validation
# folds used by all four penalized regression methods.
#
# WHY THIS FUNCTION EXISTS
# If each method drew different folds, method comparisons would include
# avoidable fold-assignment noise. A balanced, seeded fold vector makes the
# fits reproducible and directly comparable.
#
# INPUTS
# n_train: Number of training observations.
# nfolds: Number of cross-validation folds.
# seed: Seed controlling the random permutation of balanced fold labels.
#
# RETURN VALUE
# An integer vector of length n_train containing fold labels 1, ..., nfolds.
make_cv_foldid_testset <- function(n_train, nfolds, seed) {
  # Repeat labels as evenly as possible, then randomly permute them.
  set.seed(as.integer(seed))
  return(sample(
    rep(seq_len(as.integer(nfolds)), length.out = as.integer(n_train)),
    size = as.integer(n_train),
    replace = FALSE
  ))
}

# PURPOSE
# Simulate an n-by-p Gaussian predictor matrix with independent equicorrelated
# blocks of size b.
#
# WHY THIS FUNCTION EXISTS
# The test-set pilot must use the same covariance structure as the other
# simulation components while remaining standalone. The common-factor formula
# produces unit marginal variances and within-block correlation rho without
# constructing or decomposing a large p-by-p covariance matrix.
#
# INPUTS
# n: Number of observations.
# p: Number of predictors.
# b: Predictors per block.
# rho: Common correlation between two predictors in the same block.
#
# RETURN VALUE
# A numeric matrix with n rows and p columns. Different blocks are independent.
simulate_X_blockcor_testset <- function(n, p, b, rho) {
  n <- as.integer(n)
  p <- as.integer(p)
  b <- as.integer(b)
  G <- as.integer(p / b)
  # sqrt(rho) multiplies the shared block factor, sqrt(1-rho) multiplies the
  # variable-specific component. Their variances add to one.
  X <- matrix(NA_real_, nrow = n, ncol = p)
  sr <- sqrt(rho)
  si <- sqrt(1 - rho)

  # Fill one block at a time, avoiding a large covariance matrix in memory.
  for (g in seq_len(G)) {
    cols <- ((g - 1L) * b + 1L):(g * b)
    z <- rnorm(n)
    u <- matrix(rnorm(n * b), nrow = n, ncol = b)
    X[, cols] <- sr * z + si * u
  }
  return(X)
}

# PURPOSE
# Evaluate x' Sigma x for the block-equicorrelated covariance matrix without
# explicitly creating Sigma.
#
# WHY THIS FUNCTION EXISTS
# This quadratic form is needed both to scale beta to the requested population
# SNR and to calculate exact prediction-error variance. The block formula is
# faster and substantially smaller than allocating a dense p-by-p matrix,
# especially when p = 5,000.
#
# INPUTS
# x: Numeric coefficient or coefficient-error vector.
# b: Block size.
# rho: Within-block equicorrelation.
#
# RETURN VALUE
# The scalar population variance x' Sigma x.
quadform_block_equicorr_testset <- function(x, b, rho) {
  x <- as.numeric(x)
  b <- as.integer(b)
  G <- length(x) / b
  out <- 0
  # For one equicorrelation block, x' Sigma x equals
  # (1-rho) * sum(x_j^2) + rho * (sum x_j)^2.
  for (g in seq_len(G)) {
    idx <- ((g - 1L) * b + 1L):(g * b)
    xb <- x[idx]
    out <- out + (1 - rho) * sum(xb^2) + rho * sum(xb)^2
  }
  return(as.numeric(out))
}

# PURPOSE
# Generate the sparse coefficient vector for one setting and scale it to the
# requested population signal-to-noise ratio.
#
# WHY THIS FUNCTION EXISTS
# Support placement, coefficient magnitudes, signs, and SNR scaling define the
# simulated signal. Keeping them in one function ensures that training and all
# later test samples use exactly the same true beta.
#
# INPUTS
# setting_row: One row from make_testset_settings().
# seed: Seed controlling block signs and raw coefficient magnitudes.
# design: Fixed data-generating constants from make_testset_pilot_design().
#
# RETURN VALUE
# A list containing the unscaled and scaled beta vectors, their scaling
# quantities, active block/variable indices, block signs, and the minimum
# absolute scaled active coefficient retained in the saved exact-arm output.
make_beta_testset <- function(setting_row, seed,
                              design = make_testset_pilot_design()) {
  # Step 1: extract the setting and initialize all coefficients at zero.
  set.seed(as.integer(seed))

  p <- as.integer(setting_row$p)
  b <- as.integer(setting_row$b)
  B_active <- as.integer(setting_row$B_active)
  k <- as.integer(setting_row$k)
  rho <- as.numeric(setting_row$rho)
  snr <- as.numeric(setting_row$snr)

  beta0 <- rep(0, p)

  # Active blocks are deterministic: 1, ..., B_active.
  active_blocks <- seq_len(B_active)
  active_variables <- integer(0)
  block_signs <- integer(B_active)

  # Step 2: populate the deterministic active support. Randomness enters only
  # through one common sign per active block and heterogeneous magnitudes.
  for (ii in seq_along(active_blocks)) {
    g <- active_blocks[ii]
    block_cols <- ((g - 1L) * b + 1L):(g * b)

    # The first k variables in each active block are active.
    active_in_block <- block_cols[seq_len(k)]

    # Draw one common random sign for each active block.
    block_sign <- sample(c(-1L, 1L), size = 1L)
    block_signs[ii] <- block_sign

    # Draw heterogeneous positive magnitudes from U(0.5, 1.5).
    magnitudes <- runif(
      n = k,
      min = design$beta_magnitude_min,
      max = design$beta_magnitude_max
    )

    beta0[active_in_block] <- block_sign * magnitudes
    active_variables <- c(active_variables, active_in_block)
  }

  active_variables <- as.integer(active_variables)

  # Step 3: compute the raw population signal variance beta0' Sigma beta0.
  V0 <- quadform_block_equicorr_testset(beta0, b = b, rho = rho)
  if (!is.finite(V0) || V0 <= 0) {
    stop("The unscaled signal variance is not positive and finite.")
  }

  # Step 4: rescale all non-zero coefficients by one common factor so that
  # Var(X beta) / sigma^2 equals the setting's target SNR exactly.
  beta_scale <- sqrt(snr * design$sigma2 / V0)
  beta <- beta_scale * beta0

  if (any(!is.finite(beta[active_variables]))) {
    stop("A scaled active coefficient is not finite.")
  }

  return(list(
    beta0 = beta0,
    beta = beta,
    V0 = V0,
    beta_scale = beta_scale,
    active_blocks = active_blocks,
    active_variables = active_variables,
    block_signs = block_signs,
    min_abs_beta_active = min(abs(beta[active_variables]))
  ))
}

# PURPOSE
# Simulate one training data set from the Gaussian linear model for a fixed
# setting and previously generated true beta.
#
# WHY THIS FUNCTION EXISTS
# The fitted models must share the same training sample, and the test-set arm
# must later hold those fits fixed. Returning X and y together makes that
# conditioning step explicit.
#
# INPUTS
# setting_row: One setting from make_testset_settings().
# beta_obj: Output of make_beta_testset(). Its scaled beta is used here.
# seed: Seed controlling predictors and Gaussian noise.
# design: Provides sigma2.
#
# RETURN VALUE
# A list containing X_train and y_train from y = X beta + epsilon.
simulate_training_testset <- function(setting_row, beta_obj, seed,
                                      design = make_testset_pilot_design()) {
  # Generate predictors first, then independent N(0, sigma2) errors.
  set.seed(as.integer(seed))
  X <- simulate_X_blockcor_testset(
    n = setting_row$n,
    p = setting_row$p,
    b = setting_row$b,
    rho = setting_row$rho
  )
  eps <- rnorm(
    n = as.integer(setting_row$n),
    mean = 0,
    sd = sqrt(design$sigma2)
  )
  y <- as.numeric(X %*% beta_obj$beta + eps)
  return(list(X_train = X, y_train = y))
}

# PURPOSE
# Fit one cross-validated Gaussian glmnet model and extract the coefficients
# selected at lambda.min.
#
# WHY THIS FUNCTION EXISTS
# Ridge, LASSO, and both Elastic Net variants differ only in alpha. A single
# implementation guarantees identical folds, standardization, intercept
# treatment, lambda selection, and timing across methods.
#
# INPUTS
# X_train, y_train: Shared training predictors and response.
# alpha: Mixing value (0 Ridge, 1/3 EN, 0.95 near-LASSO EN, 1 LASSO).
# foldid: Shared fold allocation from make_cv_foldid_testset().
# design: glmnet configuration and selected prediction lambda.
#
# RETURN VALUE
# A list with the cv.glmnet object, estimated intercept and coefficients,
# lambda.min, lambda.1se, and elapsed fit time in seconds.
fit_one_glmnet_testset <- function(X_train, y_train, alpha, foldid,
                                   design = make_testset_pilot_design()) {
  # Time only cross-validation and fitting, coefficient extraction follows.
  t0 <- proc.time()[["elapsed"]]
  cv_fit <- glmnet::cv.glmnet(
    x = X_train,
    y = y_train,
    family = "gaussian",
    alpha = alpha,
    type.measure = design$glmnet_type_measure,
    nfolds = as.integer(design$glmnet_nfolds),
    foldid = foldid,
    standardize = design$glmnet_standardize,
    intercept = design$glmnet_intercept,
    parallel = FALSE,
    keep = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  # Extract the exact coefficient vector used for both analytical moments and
  # empirical test predictions.
  cc <- as.matrix(stats::coef(cv_fit, s = design$glmnet_prediction_s))
  intercept_hat <- as.numeric(cc[1L, 1L])
  beta_hat <- as.numeric(cc[-1L, 1L])

  return(list(
    cv_fit = cv_fit,
    intercept_hat = intercept_hat,
    beta_hat = beta_hat,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    fit_time_sec = as.numeric(elapsed)
  ))
}

# PURPOSE
# Generate one truth/training sample and fit all four penalized regression
# methods that will remain fixed during test-set evaluation.
#
# WHY THIS FUNCTION EXISTS
# The pilot isolates finite-test-set variation. Consequently beta, training
# data, folds, and fitted models must be created once and then reused across
# candidate n_test values and repeated test samples.
#
# INPUTS
# setting_row: One row from the test-set setting table.
# base_seed: Starting seed, fixed offsets separate beta, training data, and
#   cross-validation folds reproducibly.
# design: Test-set pilot design and model settings.
#
# RETURN VALUE
# A list containing beta_obj, training data, shared fold labels, and named fits
# for ridge, lasso, en, and en095.
fit_fixed_models_testset <- function(setting_row, base_seed,
                                     design = make_testset_pilot_design()) {
  # Step 1: check the glmnet dependency before simulation.
  check_required_packages_testset()

  # Step 2: create the truth and one training sample using separate seed streams.
  beta_obj <- make_beta_testset(
    setting_row = setting_row,
    seed = base_seed + 1L,
    design = design
  )
  dat <- simulate_training_testset(
    setting_row = setting_row,
    beta_obj = beta_obj,
    seed = base_seed + 2L,
    design = design
  )
  foldid <- make_cv_foldid_testset(
    n_train = nrow(dat$X_train),
    nfolds = design$glmnet_nfolds,
    seed = base_seed + 3L
  )

  # Step 3: use the same data and folds for every alpha value. This keeps
  # differences attributable to the methods rather than resampling.
  fits <- list(
    ridge = fit_one_glmnet_testset(
      dat$X_train, dat$y_train, alpha = 0, foldid = foldid, design = design
    ),
    lasso = fit_one_glmnet_testset(
      dat$X_train, dat$y_train, alpha = 1, foldid = foldid, design = design
    ),
    en = fit_one_glmnet_testset(
      dat$X_train, dat$y_train, alpha = design$en_alpha,
      foldid = foldid, design = design
    ),
    en095 = fit_one_glmnet_testset(
      dat$X_train, dat$y_train, alpha = design$en_near_lasso_alpha,
      foldid = foldid, design = design
    )
  )

  return(list(beta_obj = beta_obj, training = dat, foldid = foldid, fits = fits))
}

# PURPOSE
# Derive exact test-MSE and test-MAE moments for a fixed fitted model at every
# candidate test-set size.
#
# WHY THIS FUNCTION EXISTS
# Once beta_hat and the intercept are fixed, Gaussian predictors and noise make
# the prediction error Gaussian. Its analytical moments isolate finite-test-
# set variation without repeatedly simulating the largest test matrices.
#
# INPUTS
# beta_true: True scaled coefficient vector.
# fit_obj: One object returned by fit_one_glmnet_testset().
# b, rho: Block-covariance parameters.
# n_test_values: Candidate independent test-set sizes.
# sigma2: Response-noise variance.
#
# RETURN VALUE
# One data-frame row per n_test with the conditional error distribution,
# expected MSE/MAE, their finite-sample variances and SDs, and coefficients of
# variation.
#
# DERIVATION FOR A FIXED FITTED MODEL
# Prediction Error:
#   e = Y - Yhat = X'(beta - beta_hat) + eps - intercept_hat.
# Because both X and eps are Gaussian, the exact distribution is
#   e ~ N(mu_e, var_e)
# with
#   mu_e  = -intercept_hat,
#   var_e = sigma2 + (beta-beta_hat)' Sigma (beta-beta_hat).
# Therefore,
#   E(e^2)    = var_e + mu_e^2,
#   Var(e^2)  = 2 var_e^2 + 4 mu_e^2 var_e,
#   Var(MSE_n) = Var(e^2)/n.
conditional_test_error_moments_testset <- function(beta_true, fit_obj, b, rho,
                                                   n_test_values,
                                                   sigma2 = 1) {
  # Step 1: express prediction error through coefficient error and intercept.
  d <- beta_true - fit_obj$beta_hat
  mu_e <- -fit_obj$intercept_hat
  var_signal_error <- quadform_block_equicorr_testset(d, b = b, rho = rho)
  var_e <- sigma2 + var_signal_error

  # Step 2: use normal-distribution moments for squared prediction error.
  expected_mse <- var_e + mu_e^2
  var_sq_error <- 2 * var_e^2 + 4 * mu_e^2 * var_e

  # Step 3: absolute Gaussian errors follow a folded-normal distribution. Its
  # first two moments give the MAE expectation and variance.
  sd_e <- sqrt(var_e)
  z <- if (sd_e > 0) mu_e / sd_e else 0
  expected_mae <-
    sd_e * sqrt(2 / pi) * exp(-0.5 * z^2) +
    mu_e * (2 * stats::pnorm(z) - 1)
  var_abs_error <- expected_mse - expected_mae^2

  # Step 4: averaging n_test independent losses divides each single-loss
  # variance by n_test, expectations remain unchanged.
  return(data.frame(
    n_test = as.integer(n_test_values),
    residual_mean = mu_e,
    residual_variance = var_e,
    expected_mse = expected_mse,
    var_single_squared_error = var_sq_error,
    var_test_mse = var_sq_error / n_test_values,
    sd_test_mse = sqrt(var_sq_error / n_test_values),
    cv_test_mse = sqrt(var_sq_error / n_test_values) / expected_mse,
    expected_mae = expected_mae,
    var_single_absolute_error = var_abs_error,
    var_test_mae = var_abs_error / n_test_values,
    sd_test_mae = sqrt(var_abs_error / n_test_values),
    cv_test_mae = sqrt(var_abs_error / n_test_values) / expected_mae,
    stringsAsFactors = FALSE
  ))
}

# PURPOSE
# Fit fixed models across the test-set pilot settings and calculate exact
# finite-test-set MSE/MAE moments for every candidate n_test.
#
# WHY THIS FUNCTION EXISTS
# Repeatedly simulating test sets is unnecessary for the primary stability
# calculation: conditional on a fixed fit, the Gaussian error moments are
# known exactly. This arm therefore covers all structural settings without
# allocating a large test matrix for each of them.
#
# INPUTS
# settings: Setting table returned by make_testset_settings().
# design: Complete test-set pilot constants.
# R_base: Number of independent truth/training/fitted-model realizations per
#   setting, the documented pilot uses one.
# progress_every: Print progress after this many settings.
#
# RETURN VALUE
# A data frame with one row for every setting, base repetition, model, and
# candidate n_test. It contains exact conditional error moments plus setting,
# fit-time, lambda, and signal annotations.
run_exact_testset_variance_arm <- function(
    settings = make_testset_settings(),
    design = make_testset_pilot_design(),
    R_base = design$R_base_exact,
    progress_every = 5L) {

  # Step 1: verify glmnet and allocate a result list. Each fitted setting yields
  # four method-specific data frames, each containing all n_test candidates.
  check_required_packages_testset()
  if (length(R_base) != 1L || !is.finite(R_base) || R_base != floor(R_base) ||
    R_base < 1 || R_base > 9999) {
    stop("R_base must be one integer from 1 to 9999.")
  }
  R_base <- as.integer(R_base)
  out <- vector("list", nrow(settings) * R_base * 4L)
  pos <- 0L

  # Step 2: work through settings and optional independent training realizations.
  for (i in seq_len(nrow(settings))) {
    setting_row <- settings[i, , drop = FALSE]

    for (base_rep in seq_len(R_base)) {
      base_seed <- make_task_seeds(
        setting_row$setting_id, base_rep, design$seed_base_exact
      )

      # Create beta, training data, common folds, and the four fixed fits once.
      fixed <- fit_fixed_models_testset(
        setting_row = setting_row,
        base_seed = base_seed,
        design = design
      )

      # Step 3: derive moments at all candidate sizes separately for each fit.
      for (method in names(fixed$fits)) {
        mm <- conditional_test_error_moments_testset(
          beta_true = fixed$beta_obj$beta,
          fit_obj = fixed$fits[[method]],
          b = setting_row$b,
          rho = setting_row$rho,
          n_test_values = design$n_test_values,
          sigma2 = design$sigma2
        )

        # Attach identifiers and study factors so the row block is fully
        # interpretable after it is combined with every other setting.
        pos <- pos + 1L
        mm$setting_id <- setting_row$setting_id
        mm$base_rep <- base_rep
        mm$method <- method
        mm$n <- setting_row$n
        mm$p_over_n <- setting_row$p_over_n
        mm$p <- setting_row$p
        mm$rho <- setting_row$rho
        mm$a <- setting_row$a
        mm$k <- setting_row$k
        mm$snr <- setting_row$snr
        mm$B_active <- setting_row$B_active
        mm$q_active <- setting_row$q_active
        mm$min_abs_beta_active <- fixed$beta_obj$min_abs_beta_active
        mm$fit_time_sec <- fixed$fits[[method]]$fit_time_sec
        mm$lambda_min <- fixed$fits[[method]]$lambda_min
        out[[pos]] <- mm
      }
    }

    # Progress messages are informational only and do not affect computation.
    if (i %% progress_every == 0L || i == nrow(settings)) {
      message("Exact test-set variance arm: setting ", i, " / ", nrow(settings))
    }
  }

  # Step 4: combine all method-specific blocks into one long-format table.
  return(do.call(rbind, out[seq_len(pos)]))
}

# PURPOSE
# Aggregate exact conditional test-error moments into a compact comparison by
# predictor dimension, method, and candidate test-set size.
#
# WHY THIS FUNCTION EXISTS
# The setting-level exact output retains the standard deviation and coefficient
# of variation for every fitted setting. The design decision depends on the
# typical and worst finite-test variation, so the summary reports medians, 95th
# percentiles, and maxima.
#
# INPUT
# exact_results: Output from run_exact_testset_variance_arm().
#
# RETURN VALUE
# A data frame grouped by p, method, and n_test with MSE and MAE stability
# summaries. Rows are sorted for direct comparison across test-set sizes.
summarise_exact_testset_variance <- function(exact_results) {
  # Step 1: create one grouping key for each dimension-method-size combination.
  # Base-R split/lapply keeps the aggregation explicit and dependency-free.
  keys <- interaction(
    exact_results$p,
    exact_results$method,
    exact_results$n_test,
    drop = TRUE,
    lex.order = TRUE
  )
  spl <- split(exact_results, keys)

  # Step 2: summarize conditional results across all fitted settings in a group.
  rows <- lapply(spl, function(d) {
    return(data.frame(
      p = d$p[1],
      method = d$method[1],
      n_test = d$n_test[1],
      n_fit_settings = nrow(d),
      median_sd_test_mse = median(d$sd_test_mse),
      q95_sd_test_mse = as.numeric(quantile(d$sd_test_mse, 0.95)),
      max_sd_test_mse = max(d$sd_test_mse),
      median_cv_test_mse = median(d$cv_test_mse),
      q95_cv_test_mse = as.numeric(quantile(d$cv_test_mse, 0.95)),
      max_cv_test_mse = max(d$cv_test_mse),
      median_sd_test_mae = median(d$sd_test_mae),
      median_cv_test_mae = median(d$cv_test_mae),
      stringsAsFactors = FALSE
    ))
  })

  # Step 3: bind and sort the small group summaries into one readable table.
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out[order(out$p, out$method, out$n_test), , drop = FALSE])
}

# PURPOSE
# Select one structurally demanding setting for each of the two predictor
# dimensions used in the empirical runtime arm.
#
# WHY THIS FUNCTION EXISTS
# Matrix generation and prediction time are driven mainly by p and n_test.
# Benchmarking every structural factor combination would largely repeat the
# same size calculation. The runtime arm therefore uses one high-correlation,
# dense-within-active-block reference setting at each dimension.
#
# INPUT
# settings: Full table from make_testset_settings().
#
# RETURN VALUE
# The two retained setting rows (p = 40 and p = 5,000).
select_runtime_benchmark_settings <- function(
    settings = make_testset_settings()) {
  # Use one structurally demanding setting per dimensional extreme. Runtime
  # depends mainly on p and n_test. The values of rho, k, and SNR do not change
  # the matrix dimensions.
  # The common benchmark uses rho = .95, a = .75, k = 10, and SNR = 2.
  keep <- settings[
    settings$rho == 0.95 &
      settings$a == 0.75 &
      settings$k == 10L &
      settings$snr == 2,
    , drop = FALSE
  ]
  return(keep)
}

# PURPOSE
# Evaluate four already fitted models on one nested sequence of independent
# test samples and measure test-data/prediction time.
#
# WHY THIS FUNCTION EXISTS
# A single simulated sequence makes estimates at n_test = 2,500, 5,000,
# 10,000, and 20,000 nested and directly comparable. Chunking prevents the
# p = 5,000 arm from allocating the largest possible test matrix at once while
# leaving the statistical sample unchanged.
#
# INPUTS
# setting_row: One selected runtime benchmark setting.
# fixed: Output from fit_fixed_models_testset(), models are not refitted here.
# test_seed: Seed for the independent test sequence.
# design: Candidate sizes, chunk size, noise variance, and prediction lambda.
#
# RETURN VALUE
# A data frame with one row per method and n_test checkpoint. It contains
# cumulative MSE, MAE, predictor/response generation time, method-specific
# prediction time, metric time, and total elapsed time.
run_nested_test_evaluation_testset <- function(setting_row, fixed,
                                               test_seed,
                                               design = make_testset_pilot_design()) {
  # Step 1: initialize cumulative loss sums and timing counters. Cumulative
  # sums let each larger checkpoint extend, rather than replace, the smaller one.
  checkpoints <- sort(as.integer(design$n_test_values))
  max_n <- max(checkpoints)
  chunk_size <- as.integer(design$runtime_chunk_size)
  methods <- names(fixed$fits)

  sse <- setNames(rep(0, length(methods)), methods)
  sae <- setNames(rep(0, length(methods)), methods)
  time_pred <- setNames(rep(0, length(methods)), methods)
  time_X <- 0
  time_y <- 0
  time_metric <- 0
  n_done <- 0L
  result_rows <- list()
  result_pos <- 0L

  # Step 2: one seed defines the complete nested test sequence.
  set.seed(as.integer(test_seed))
  total_start <- proc.time()[["elapsed"]]

  # Step 3: generate, predict, and discard one manageable chunk at a time.
  while (n_done < max_n) {
    n_chunk <- min(chunk_size, max_n - n_done)

    tx <- proc.time()[["elapsed"]]
    X <- simulate_X_blockcor_testset(
      n = n_chunk,
      p = setting_row$p,
      b = setting_row$b,
      rho = setting_row$rho
    )
    time_X <- time_X + (proc.time()[["elapsed"]] - tx)

    ty <- proc.time()[["elapsed"]]
    eps <- rnorm(n_chunk, mean = 0, sd = sqrt(design$sigma2))
    y <- as.numeric(X %*% fixed$beta_obj$beta + eps)
    time_y <- time_y + (proc.time()[["elapsed"]] - ty)

    # Apply every fixed model to exactly the same new observations.
    for (method in methods) {
      tp <- proc.time()[["elapsed"]]
      pred <- as.numeric(stats::predict(
        fixed$fits[[method]]$cv_fit,
        newx = X,
        s = design$glmnet_prediction_s,
        type = "response"
      ))
      time_pred[method] <- time_pred[method] +
        (proc.time()[["elapsed"]] - tp)

      tm <- proc.time()[["elapsed"]]
      err <- y - pred
      sse[method] <- sse[method] + sum(err^2)
      sae[method] <- sae[method] + sum(abs(err))
      time_metric <- time_metric + (proc.time()[["elapsed"]] - tm)
    }

    n_done <- n_done + n_chunk

    # Step 4: whenever a requested sample size is reached, convert cumulative
    # sums into MSE/MAE and record cumulative timings.
    if (n_done %in% checkpoints) {
      total_elapsed <- proc.time()[["elapsed"]] - total_start
      for (method in methods) {
        result_pos <- result_pos + 1L
        result_rows[[result_pos]] <- data.frame(
          n_test = n_done,
          method = method,
          mse = sse[method] / n_done,
          mae = sae[method] / n_done,
          cumulative_time_X_sec = time_X,
          cumulative_time_y_sec = time_y,
          cumulative_time_prediction_sec = time_pred[method],
          cumulative_time_metrics_all_methods_sec = time_metric,
          cumulative_time_total_sec = total_elapsed,
          stringsAsFactors = FALSE
        )
      }
    }

    # The next chunk does not require the previous matrix. Releasing it keeps
    # peak memory governed by runtime_chunk_size rather than max(n_test).
    rm(X, eps, y)
    gc(FALSE)
  }

  return(do.call(rbind, result_rows))
}

# PURPOSE
# Run the empirical test-generation and prediction benchmark at p = 40 and
# p = 5,000 while holding one trained set of models fixed at each dimension.
#
# WHY THIS FUNCTION EXISTS
# The analytical arm answers the statistical stability question but does not
# measure the practical cost of creating large predictor matrices and making
# predictions. This complementary arm records that cost for the two matrix
# extremes without timing model fitting repeatedly.
#
# INPUTS
# settings: Usually the two rows from select_runtime_benchmark_settings().
# design: Test-set pilot constants.
# R_eval: Number of independent test sequences evaluated for each dimension.
#
# RETURN VALUE
# A row-bound data frame containing checkpoint metrics and cumulative timings,
# annotated with repetition, setting, and factor values.
run_runtime_testset_arm <- function(
    settings = select_runtime_benchmark_settings(),
    design = make_testset_pilot_design(),
    R_eval = design$R_runtime_eval) {

  # Step 1: check glmnet and initialize the result collector.
  check_required_packages_testset()
  if (length(R_eval) != 1L || !is.finite(R_eval) || R_eval != floor(R_eval) ||
    R_eval < 1 || R_eval > 9999) {
    stop("R_eval must be one integer from 1 to 9999.")
  }
  rows <- list()
  pos <- 0L

  # Step 2: fit each dimensional benchmark once.
  for (i in seq_len(nrow(settings))) {
    setting_row <- settings[i, , drop = FALSE]
    base_seed <- make_task_seeds(
      setting_row$setting_id, 1L, design$seed_base_runtime
    )

    # Fit each p benchmark once, only the test data vary afterward.
    fixed <- fit_fixed_models_testset(
      setting_row = setting_row,
      base_seed = base_seed,
      design = design
    )

    # Step 3: reuse the fixed models while drawing independent nested test
    # sequences. This isolates test-evaluation time and variability.
    for (rr in seq_len(as.integer(R_eval))) {
      # Independent evaluation families keep repeated test samples separate
      # from the fixed truth, training data, and shared folds.
      test_seed <- make_task_seeds(
        setting_row$setting_id, rr, design$seed_base_evaluation
      ) + 2L
      z <- run_nested_test_evaluation_testset(
        setting_row = setting_row,
        fixed = fixed,
        test_seed = test_seed,
        design = design
      )
      z$runtime_rep <- rr
      z$setting_id <- setting_row$setting_id
      z$n <- setting_row$n
      z$p <- setting_row$p
      z$rho <- setting_row$rho
      z$a <- setting_row$a
      z$k <- setting_row$k
      z$snr <- setting_row$snr
      pos <- pos + 1L
      rows[[pos]] <- z
    }
  }

  return(do.call(rbind, rows))
}

# PURPOSE
# Summarize the empirical runtime arm by dimension, method, and test-set size.
#
# WHY THIS FUNCTION EXISTS
# Individual evaluations show run-to-run variation, whereas a design decision
# needs a compact comparison of average time and empirical stability. The
# summary therefore reports standard deviations and coefficients of variation,
# not individual or mean prediction-error values.
#
# INPUT
# runtime_results: Output of run_runtime_testset_arm().
#
# RETURN VALUE
# A sorted data frame with the evaluation count, empirical standard deviations
# and coefficients of variation for both prediction-error measures, and mean
# cumulative generation, prediction, and total times.
summarise_runtime_testset_arm <- function(runtime_results) {
  # Step 1: split repetitions into comparable dimension-method-size groups.
  keys <- interaction(
    runtime_results$p,
    runtime_results$method,
    runtime_results$n_test,
    drop = TRUE,
    lex.order = TRUE
  )
  spl <- split(runtime_results, keys)
  # Step 2: aggregate repeated nested evaluations within each group.
  rows <- lapply(spl, function(d) {
    mean_mse <- mean(d$mse)
    mean_mae <- mean(d$mae)
    sd_mse <- if (nrow(d) > 1) sd(d$mse) else NA_real_
    sd_mae <- if (nrow(d) > 1) sd(d$mae) else NA_real_

    return(data.frame(
      p = d$p[1],
      method = d$method[1],
      n_test = d$n_test[1],
      R_eval = nrow(d),
      sd_mse_empirical = sd_mse,
      cv_mse_empirical = sd_mse / mean_mse,
      sd_mae_empirical = sd_mae,
      cv_mae_empirical = sd_mae / mean_mae,
      mean_time_X_sec = mean(d$cumulative_time_X_sec),
      mean_time_y_sec = mean(d$cumulative_time_y_sec),
      mean_time_prediction_sec = mean(d$cumulative_time_prediction_sec),
      mean_time_total_sec = mean(d$cumulative_time_total_sec),
      stringsAsFactors = FALSE
    ))
  })
  # Step 3: bind and sort the summaries for a readable size comparison.
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out[order(out$p, out$method, out$n_test), , drop = FALSE])
}

# PURPOSE
# Show the theoretical 1/n_test variance and 1/sqrt(n_test)
# standard-deviation scaling relative to the selected reference size of 5,000.
#
# WHY THIS FUNCTION EXISTS
# The table provides a model-independent interpretation of what doubling or
# halving a test set can achieve. The final column gives the MSE coefficient of
# variation for zero-mean Gaussian prediction error as a simple benchmark.
#
# INPUT
# n_test_values: Candidate independent test-set sizes.
#
# RETURN VALUE
# A data frame with variance and standard-deviation ratios relative to
# n_test = 5,000 and the Gaussian-error reference coefficient of variation.
make_theoretical_n_test_scaling_table <- function(
    n_test_values = make_testset_pilot_design()$n_test_values) {
  n_test_values <- sort(as.integer(n_test_values))
  ref <- 5000
  return(data.frame(
    n_test = n_test_values,
    var_ratio_vs_5000 = ref / n_test_values,
    sd_ratio_vs_5000 = sqrt(ref / n_test_values),
    cv_mse_if_zero_mean_gaussian_error = sqrt(2 / n_test_values),
    stringsAsFactors = FALSE
  ))
}

# PURPOSE
# Prepare the four public result tables for the test-set-size pilot.
#
# WHY THIS FUNCTION EXISTS
# The analytical and runtime calculations need prediction-error means and
# individual evaluation values internally. Those values are intermediate
# quantities, not outputs needed to justify the selected test-set size. Keeping
# the public column selection in one function ensures that comma-separated
# values files, the summary R data object, and the object returned by the
# launcher follow the same rule. Public results retain standard deviations,
# coefficients of variation, design information, tuning values, and timings.
#
# INPUTS
# exact_results: Detailed analytical results created internally by
#   run_exact_testset_variance_arm().
# exact_summary: Aggregated analytical results created by
#   summarise_exact_testset_variance().
# runtime_results: Detailed empirical timing results created internally by
#   run_runtime_testset_arm().
# runtime_summary: Aggregated timing and variability results created by
#   summarise_runtime_testset_arm().
#
# RETURN VALUE
# A named list of four data frames. No table contains an individual or mean
# prediction mean squared error or prediction mean absolute error.
make_testset_public_tables <- function(exact_results,
                                       exact_summary,
                                       runtime_results,
                                       runtime_summary) {
  exact_columns <- c(
    "n_test", "sd_test_mse", "cv_test_mse", "sd_test_mae",
    "cv_test_mae", "setting_id", "base_rep", "method", "n",
    "p_over_n", "p", "rho", "a", "k", "snr", "B_active",
    "q_active", "min_abs_beta_active", "fit_time_sec", "lambda_min"
  )

  runtime_columns <- setdiff(names(runtime_results), c("mse", "mae"))

  return(list(
    exact_results = exact_results[, exact_columns, drop = FALSE],
    exact_summary = exact_summary,
    runtime_results = runtime_results[, runtime_columns, drop = FALSE],
    runtime_summary = runtime_summary
  ))
}

# PURPOSE
# Save the public numerical products of the test-set-size pilot in readable
# comma-separated values tables and one summary R data object.
#
# WHY THIS FUNCTION EXISTS
# Comma-separated values files support direct inspection and reuse outside R,
# while the R data file preserves types, the design, and the same public
# tables in one object. The inputs have already passed through
# make_testset_public_tables(), so individual and mean prediction errors cannot
# enter either format. Timestamped names distinguish separate executions.
#
# INPUTS
# exact_results, exact_summary: Detailed and aggregated analytical-arm tables.
# runtime_results, runtime_summary: Detailed and aggregated runtime-arm tables.
# output_dir: Destination directory.
# run_stamp: Timestamp shared by every file from this run. The orchestration
#   function creates it once and also returns it for the session-info filename.
#
# RETURN VALUE
# A named character vector containing the paths of all six written files.
save_testset_pilot_outputs <- function(exact_results,
                                       exact_summary,
                                       runtime_results,
                                       runtime_summary,
                                       output_dir = file.path(
                                         "outputs", "testset_size_pilot"
                                       ),
                                       run_stamp = format(
                                         Sys.time(), "%Y-%m-%d_%H-%M-%S"
)) {
  # Step 1: create the destination and use the timestamp supplied for this run.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- run_stamp

  f_exact <- file.path(output_dir, paste0("testset_exact_results_", stamp, ".csv"))
  f_exact_sum <- file.path(output_dir, paste0("testset_exact_summary_", stamp, ".csv"))
  f_runtime <- file.path(output_dir, paste0("testset_runtime_results_", stamp, ".csv"))
  f_runtime_sum <- file.path(output_dir, paste0("testset_runtime_summary_", stamp, ".csv"))
  f_scaling <- file.path(output_dir, paste0("testset_theoretical_scaling_", stamp, ".csv"))
  f_rds <- file.path(output_dir, paste0("testset_size_pilot_", stamp, ".rds"))

  # Step 2: write flat tables that can be inspected without loading R objects.
  utils::write.csv(exact_results, f_exact, row.names = FALSE)
  utils::write.csv(exact_summary, f_exact_sum, row.names = FALSE)
  utils::write.csv(runtime_results, f_runtime, row.names = FALSE)
  utils::write.csv(runtime_summary, f_runtime_sum, row.names = FALSE)
  utils::write.csv(make_theoretical_n_test_scaling_table(), f_scaling, row.names = FALSE)

  saved_files <- c(
    exact_results = f_exact,
    exact_summary = f_exact_sum,
    runtime_results = f_runtime,
    runtime_summary = f_runtime_sum,
    theoretical_scaling = f_scaling,
    rds = f_rds
  )

  # Step 3: store the same public summaries in one R object. This object is
  # convenient for R users but does not add hidden repetition-level results.
  # The timestamp and file map remain available if the file is moved.
  saveRDS(
    list(
      run_stamp = stamp,
      design = make_testset_pilot_design(),
      exact_results = exact_results,
      exact_summary = exact_summary,
      runtime_results = runtime_results,
      runtime_summary = runtime_summary,
      theoretical_scaling = make_theoretical_n_test_scaling_table(),
      saved_files = saved_files,
      created_at = Sys.time()
    ),
    f_rds
  )

  return(saved_files)
}

# PURPOSE
# Orchestrate the complete test-set-size pilot from setting construction to
# exact moments, empirical runtime measurements, summaries, and saved files.
#
# WHY THIS FUNCTION EXISTS
# The user-facing launcher needs one clear entry point that runs the two distinct
# arms in the correct order and returns a self-contained object. This function
# does not participate in selecting the study-wide Monte Carlo repetition count.
#
# INPUTS
# output_dir: Destination for all pilot products.
# R_base_exact: Number of fixed training fits per setting in the analytical arm.
# R_runtime_eval: Number of independent nested test sequences per runtime
#   benchmark setting.
#
# RETURN VALUE
# A list containing the shared run_stamp, design, settings, public variability
# and timing tables from both arms, the theoretical scaling table, and paths of
# all saved files. Individual and mean prediction errors are not returned.
run_testset_size_pilot <- function(
    output_dir = file.path("outputs", "testset_size_pilot"),
    R_base_exact = make_testset_pilot_design()$R_base_exact,
    R_runtime_eval = make_testset_pilot_design()$R_runtime_eval) {

  # Step 1: create the fixed design and the two-dimensional-extreme setting grid.
  design <- make_testset_pilot_design()
  settings <- make_testset_settings(design)

  # Step 2: fit models and derive conditional MSE/MAE moments analytically.
  message("Test-set-size pilot: starting the exact variance arm.")
  exact_results <- run_exact_testset_variance_arm(
    settings = settings,
    design = design,
    R_base = R_base_exact
  )
  exact_summary <- summarise_exact_testset_variance(exact_results)

  # Step 3: separately benchmark chunked test generation and prediction.
  message("Test-set-size pilot: starting the separate runtime arm.")
  runtime_results <- run_runtime_testset_arm(
    settings = select_runtime_benchmark_settings(settings),
    design = design,
    R_eval = R_runtime_eval
  )
  runtime_summary <- summarise_runtime_testset_arm(runtime_results)

  # Step 4: remove individual and mean prediction-error values after all
  # standard deviations and coefficients of variation have been calculated.
  # The raw objects remain only in this function's temporary workspace and are
  # released when the function returns.
  public_tables <- make_testset_public_tables(
    exact_results = exact_results,
    exact_summary = exact_summary,
    runtime_results = runtime_results,
    runtime_summary = runtime_summary
  )

  # Step 5: create one run identifier and write all public products with that
  # timestamp. The launcher uses the same identifier for the software-session
  # record.
  run_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  saved <- save_testset_pilot_outputs(
    exact_results = public_tables$exact_results,
    exact_summary = public_tables$exact_summary,
    runtime_results = public_tables$runtime_results,
    runtime_summary = public_tables$runtime_summary,
    output_dir = output_dir,
    run_stamp = run_stamp
  )

  return(list(
    run_stamp = run_stamp,
    design = design,
    settings = settings,
    exact_results = public_tables$exact_results,
    exact_summary = public_tables$exact_summary,
    runtime_results = public_tables$runtime_results,
    runtime_summary = public_tables$runtime_summary,
    theoretical_scaling = make_theoretical_n_test_scaling_table(
      design$n_test_values
    ),
    saved_files = saved
  ))
}

###############################################################
# EXECUTION
###############################################################
# Sourcing this file only defines functions. The recommended entry point is
# run/01_run_testset_size_pilot.R. For a custom run after package installation
# and tests, call:
#
# testset_result <- run_testset_size_pilot(
#   output_dir = file.path("outputs", "testset_size_pilot"),
#   R_base_exact = 1L,
#   R_runtime_eval = 5L
# )
#
# testset_result$theoretical_scaling
# testset_result$exact_summary
# testset_result$runtime_summary
###############################################################
