###############################################################
# SHARED PILOT CORE FOR THE SIMULATION PILOTS
###############################################################
#
# PURPOSE
# This file contains the statistical and computational definitions shared by
# the precision and runtime pilots. It defines the factorial design, Gaussian
# data-generating mechanism, prediction methods, performance measures,
# deterministic random-number seeds, task tables, and serial or batchtools
# execution. The separate test-set-size pilot contains a smaller implementation
# of the same data-generating mechanism. A focused verification test compares
# that implementation with this core.
#
# WHY THIS FILE EXISTS
# The pilot scripts answer different planning questions, but they must be
# based on the same experiment. The precision pilot estimates Monte Carlo
# variability without evaluating the confirmatory Ridge comparison. The
# runtime pilot measures the computing cost of a full study repetition. This
# shared core makes those differences explicit while reusing the same
# data-generating mechanism and model code. This package contains planning
# simulations only.
#
# WHEN AND HOW TO USE IT
# Do not normally run this file on its own. It is sourced automatically by
# R/precision_pilot.R and R/runtime_pilot.R. For a new environment, use this
# order:
#   1. Execute run/00_install_packages.R once.
#   2. Execute tests/run_all_tests.R to verify the installation.
#   3. Execute the required launcher in run/, numbered 01, 02, 03, or 04.
#
# STUDY COMPONENTS DEFINED HERE
# - the complete factorial pilot design with a = {0.25, 0.75},
# - block-correlated Gaussian predictors with block size b = 10,
# - deterministic active support defined by a and k,
# - heterogeneous raw coefficient magnitudes from U(0.5, 1.5),
# - one randomly sampled common sign per active block,
# - exact population signal-to-noise ratio scaling based on the theoretical
#   covariance matrix,
# - a common independent test-set size of 5,000 observations,
# - Ridge, LASSO, Elastic Net(1/3), and Elastic Net(0.95),
# - prediction mean squared error as the primary metric, with prediction mean
#   absolute error and coefficient mean squared errors as secondary metrics,
# - the primary relative prediction-error estimand
#       theta_RL = 1 - E(MSE_Ridge) / E(MSE_Lasso)
#   and its paired influence-value Monte Carlo standard error.
#
# KEY TERMS USED THROUGHOUT THE CODE
# - design: the fixed factor levels and method settings for a pilot
# - setting: one unique combination of n, p/n, rho, a, k, and SNR
# - n: training-sample size
# - p: number of predictor variables
# - p/n: dimensionality relative to the training-sample size
# - rho: correlation between predictors within the same block
# - a: proportion of predictor blocks that contain active variables
# - k: number of active variables inside every active block
# - SNR: population signal-to-noise ratio before a data set is sampled
# - repetition: one independent realization of beta, data, and cross-validation
#   folds within a setting
# - task: one setting-repetition pair submitted to a serial or parallel worker
# - task result: the metrics and timing information returned by that task
# - setting summary: quantities obtained by aggregating repetitions within one
#   setting
#
# These levels are kept distinct throughout the code. The design contains all
# settings, each setting is repeated independently, and each completed task
# contributes one result row before setting and global summaries are formed.
# For example, n_test is the number of observations within one repetition,
# whereas R is the number of independent repetitions for every setting.
#
# NAMING OF PREDICTION-ERROR QUANTITIES
# The internal fields are named mse and mae because short general-purpose
# helper functions calculate them. Since evaluation uses an independently
# generated test sample, these quantities are prediction mean squared error
# and prediction mean absolute error. Individual values are used internally
# when a summary requires them, but they are not part of the public result
# tables. The public files contain the relevant standard deviations,
# coefficients of variation, dependence measures, and planning summaries.
# An estimand is the population-level quantity that the simulation is designed
# to estimate. Here the primary estimand compares the expected prediction mean
# squared errors of two methods through a ratio of means.
#
# FILE STRUCTURE
#   A  packages and design constants
#   B  general helper functions
#   C  factorial settings and the data-generating mechanism
#   D  shared cross-validation folds, fitted models, and performance measures
#   E  ratio-of-means estimand, influence values, and the repetition formula
#   F  reproducible task tables and one-task worker functions
#   G  detailed runtime and scheduler summaries
#   H  serial and batchtools execution, signatures, and registry locks
#   I  usage notes
#
# Scientific constants belong in make_design(). Functions receive the design
# explicitly so that the assumptions used in every calculation remain visible.
#
# METHODOLOGICAL REFERENCES USED IN THE REGISTERED REPORT
# Morris et al. (2019) describe the aims, data-generating
# mechanisms, estimands, methods, and performance-measure structure for
# simulation studies and emphasize Monte Carlo uncertainty.
#
# Siepe et al. (2024) provide a standardized framework for planning,
# preregistering, and reporting methodological simulation studies.
#
# Deng et al. (2018) describe the delta method for ratio metrics.
# That result supports the influence-value calculation used for the
# ratio-of-means estimand below.
#
# Numerical choices such as the Monte Carlo standard-error target, the minimum
# repetition count, the Random Forest sensitivity multipliers, and the
# eight-worker setup are study-specific decisions recorded in the Registered
# Report rather than values taken from these references.
###############################################################

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  RANGER_NUM_THREADS = "1"
)

###############################################################
# A. PACKAGES AND DESIGN
###############################################################

required_packages <- c(
  "batchtools", "dplyr", "tibble", "glmnet", "ranger", "snow"
)

suppressPackageStartupMessages({
  library(batchtools)
  library(dplyr)
  library(tibble)
})

# PURPOSE
# Create the one design object passed to all shared pilot functions.
#
# WHY IT EXISTS
# Keeping the scientific constants in one object prevents the precision and
# runtime pilots from silently using different factor levels, model settings,
# or test-set sizes. Only constants used by these pilots are retained here.
# Confirmatory hypothesis constants and unused planning diagnostics are
# intentionally excluded.
#
# INPUTS
# None. This package represents one fixed pilot design.
#
# RETURNS
# A named list containing factor levels, data-generating constants, model
# settings, precision targets, and the batchtools backend choice.
make_design <- function() {
  design <- list()
  design$design_name <- "ridge_lasso_gaussian_linear_pilots"

  # n and p/n jointly determine p. The remaining factors control correlation,
  # signal prevalence across blocks, signal density within active blocks, and
  # population signal strength.
  design$n_values <- c(200, 1000)
  design$p_over_n_values <- c(0.2, 5)
  design$rho_values <- c(0, 0.6, 0.8, 0.95)
  design$a_values <- c(0.25, 0.75)
  design$k_values <- c(1, 5, 10)
  design$snr_values <- c(0.5, 2)

  # Data-generating mechanism.
  design$b <- 10L
  design$sigma2 <- 1
  design$beta_magnitude_min <- 0.5
  design$beta_magnitude_max <- 1.5

  # Precision-planning target. The pilot selects R so that the estimated Monte
  # Carlo standard error does not exceed this value, subject to the stated
  # minimum. These numerical values are prespecified study choices rather than
  # values taken from Morris et al. (2019).
  design$mcse_target <- 0.002
  design$R_min <- 1000L

  # The separate test-set-size pilot supports one common size of 5,000.
  design$n_test <- 5000L
  design$n_test_policy <- "constant_5000_all_settings"

  # Random Forest is used only as a blinded precision-pilot surrogate.
  design$rf_num_trees <- 200L
  design$rf_min_node_size <- 5L

  # Penalized-regression settings passed to glmnet.
  design$glmnet_nfolds <- 5L
  design$glmnet_type_measure <- "mse"
  design$glmnet_standardize <- TRUE
  design$glmnet_intercept <- TRUE
  design$glmnet_prediction_s <- "lambda.min"
  # glmnet weights L1 by alpha and squared L2 by (1-alpha)/2.
  # alpha = 1/3 gives equal coefficients on these two penalty terms.
  design$en_alpha <- 1 / 3
  design$en_near_lasso_alpha <- 0.95

  # batchtools uses socket workers because this also works on Windows and on
  # the server configuration used for the completed pilots.
  design$batchtools_cluster_function <- "socket"

  return(design)
}

# Create the documented default once for launchers and interactive inspection.
# Functions still accept an explicit design argument, which lets the small test
# suites reduce selected dimensions without changing this package default.
design <- make_design()

###############################################################
# B. GENERAL HELPERS
###############################################################

# PURPOSE
# Produce a filesystem-safe timestamp for result filenames.
#
# WHY IT EXISTS
# Pilot outputs are expensive to regenerate. Adding a common timestamp allows
# products from separate executions to remain distinguishable.
#
# INPUTS
# None.
#
# RETURNS
# A character string such as "2026-09-22_21-30-00".
make_timestamp <- function() {
  return(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
}

# PURPOSE
# Identify the R file that is currently being executed or sourced.
#
# WHY IT EXISTS
# batchtools workers must source the shared code in a fresh R process. This
# helper supplies the absolute path without assuming a particular working
# directory or parent-folder name.
#
# INPUTS
# None. It inspects command-line arguments and active source/sys.source frames.
#
# RETURNS
# The normalized script path, or character(0) if no script path is available.
get_current_script_path <- function() {
  # Only accept the core itself. batchtools uses sys.source(), which records
  # its path as file rather than ofile. An outer test or launcher frame must
  # never be mistaken for the core when locating its sibling dependency.
  frames <- sys.frames()
  paths <- unlist(lapply(frames, function(frame) {
    return(Filter(function(path) {
      return(is.character(path) && length(path) == 1L && !is.na(path) &&
        basename(path) == "pilot_core.R")
    }, list(frame$ofile, frame$file)))
  }), use.names = FALSE)
  paths <- paths[file.exists(paths)]
  if (length(paths) > 0L) {
    return(normalizePath(tail(paths, 1L), winslash = "/", mustWork = TRUE))
  }

  # If this file is executed directly with Rscript rather than sourced, there
  # is no active source frame and the --file argument is the correct fallback.
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", cmd, fixed = TRUE)
  if (length(hit) > 0L) {
    command_file <- sub("--file=", "", cmd[hit[1]], fixed = TRUE)
    if (basename(command_file) == "pilot_core.R" && file.exists(command_file)) {
      return(normalizePath(command_file, winslash = "/", mustWork = TRUE))
    }
  }
  # Interactive loading may have no source frame. Use only explicit core paths.
  candidates <- c(file.path("R", "pilot_core.R"), "pilot_core.R")
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(normalizePath(existing[1L], winslash = "/", mustWork = TRUE))
  }
  return(character(0))
}

current_script_file <- get_current_script_path()

# Load base-R checks beside this core, also when a socket worker sources a
# copied core file from another working directory.
reproducibility_script_file <- normalizePath(
  file.path(dirname(current_script_file), "reproducibility.R"),
  winslash = "/", mustWork = TRUE
)
source(reproducibility_script_file)
check_reference_environment(packages = c(required_packages, "Matrix"))

# PURPOSE
# Calculate a mean after discarding non-finite values.
#
# WHY IT EXISTS
# This helper summarizes finite runtime measurements and returns NA_real_
# for an empty finite subset.
#
# INPUTS
# x: Numeric vector.
#
# RETURNS
# One numeric mean, or NA_real_ when no finite value is available.
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  return(if (length(x) == 0) NA_real_ else mean(x))
}

# PURPOSE
# Calculate a median after discarding non-finite values.
#
# WHY IT EXISTS
# This helper summarizes finite timing values consistently throughout
# the runtime summaries, returning NA_real_ for an empty finite subset.
#
# INPUTS
# x: Numeric vector.
#
# RETURNS
# One numeric median, or NA_real_ when no finite value is available.
safe_median <- function(x) {
  x <- x[is.finite(x)]
  return(if (length(x) == 0) NA_real_ else stats::median(x))
}

# PURPOSE
# Calculate a sample standard deviation after discarding non-finite values.
#
# WHY IT EXISTS
# At least two valid repetitions are required to estimate variability. This
# helper applies that rule consistently and returns NA_real_ for shorter subsets.
#
# INPUTS
# x: Numeric vector.
#
# RETURNS
# One sample standard deviation, or NA_real_ when fewer than two finite values
# are available.
safe_sd <- function(x) {
  x <- x[is.finite(x)]
  return(if (length(x) < 2) NA_real_ else stats::sd(x))
}

# PURPOSE
# Find the largest finite value in a numeric vector.
#
# WHY IT EXISTS
# Global precision planning takes maxima across finite setting summaries.
# This helper returns the maximum, or NA_real_ for an empty finite subset.
#
# INPUTS
# x: Numeric vector.
#
# RETURNS
# The maximum finite value, or NA_real_ when none is available.
safe_max <- function(x) {
  x <- x[is.finite(x)]
  return(if (length(x) == 0) NA_real_ else max(x))
}

###############################################################
# C. SETTINGS AND DATA-GENERATING MECHANISM
###############################################################

# PURPOSE
# Convert one combination of experimental factors into a complete setting row.
#
# WHY IT EXISTS
# Later functions should not repeatedly derive p, the number of blocks, or the
# number of active coefficients. Computing them once makes every task row
# self-contained and easy to inspect.
#
# INPUTS
# n: Training-sample size.
# p_over_n: Dimensionality ratio p/n.
# rho: Within-block predictor correlation.
# a: Proportion of blocks that contain signal.
# k: Number of active variables in every active block.
# snr: Target population signal-to-noise ratio.
# design: Shared design list returned by make_design().
#
# RETURNS
# A one-row tibble containing the original factors and the derived quantities
# p, G (total blocks), active/inactive block counts, active-variable count,
# realized block/variable sparsity, and n_test. These descriptive quantities
# are retained because they are written to the task and setting outputs.
make_setting <- function(n, p_over_n, rho, a, k, snr,
                         design = make_design()) {
  # Convert p/n to the actual number of predictors and divide those predictors
  # into equally sized blocks.
  b <- as.integer(design$b)
  p <- as.integer(round(n * p_over_n))
  G <- as.integer(p / b)

  # The first B_active blocks contain signal. Each contributes exactly k
  # active coefficients, so q_active is the total support size.
  B_active <- as.integer(round(a * G))
  q_active <- as.integer(B_active * k)

  return(tibble(
    n = as.integer(n),
    p_over_n = as.numeric(p_over_n),
    p = p,
    rho = as.numeric(rho),
    a = as.numeric(a),
    a_realized = B_active / G,
    k = as.integer(k),
    block_density = as.numeric(k / b),
    snr = as.numeric(snr),
    b = b,
    G = G,
    B_active = B_active,
    B_inactive = G - B_active,
    q_active = q_active,
    s_realized = q_active / p,
    n_test = as.integer(design$n_test)
  ))
}

# PURPOSE
# Build the complete factorial table used by the shared pilots.
#
# WHY IT EXISTS
# Precision and runtime tasks must be generated from the same ordered design.
# Building the grid in one function prevents duplicated factor definitions and
# gives every setting a stable identifier used in seeds and saved results.
#
# INPUTS
# design: Shared design list containing all factor levels.
#
# RETURNS
# One tibble with one row per factor combination and a stable setting_id.
#
# IMPLEMENTATION NOTE
# Explicit loops make the nesting order transparent. setting_id therefore has
# a deterministic meaning and is reproduced exactly whenever the design is
# rebuilt.
make_all_settings <- function(design = make_design()) {
  rows <- list()
  ii <- 1L
  for (n in design$n_values) {
    for (p_over_n in design$p_over_n_values) {
      for (rho in design$rho_values) {
        for (a in design$a_values) {
          for (k in design$k_values) {
            for (snr in design$snr_values) {
              rows[[ii]] <- make_setting(
                n, p_over_n, rho, a, k, snr, design
              )
              ii <- ii + 1L
            }
          }
        }
      }
    }
  }
  out <- bind_rows(rows)
  out$setting_id <- seq_len(nrow(out))
  out <- out %>% relocate(setting_id)
  return(out)
}

# PURPOSE
# Simulate an n x p Gaussian predictor matrix with equicorrelated blocks.
#
# WHY IT EXISTS
# Both pilots require predictors with the report's exact covariance structure.
# The common-factor construction generates that distribution directly without
# allocating and decomposing a dense p-by-p covariance matrix.
#
# STATISTICAL IDEA
# For variable j in block g,
#   X_ij = sqrt(rho) * Z_ig + sqrt(1-rho) * U_ij,
# where all latent variables are independent standard normal. Variables in the
# same block consequently have correlation rho, variables in different blocks
# are independent.
#
# INPUTS
# n: Number of observations.
# p: Number of predictors.
# b: Number of predictors per block.
# rho: Correlation shared by every pair within a block.
#
# RETURNS
# Numeric n x p matrix with columns named x1, ..., xp.
simulate_X_blockcor <- function(n, p, b, rho) {
  n <- as.integer(n)
  p <- as.integer(p)
  b <- as.integer(b)
  # At rho = 0, no common block factor is needed. Direct generation is faster
  # and produces the same independent standard-normal distribution.
  if (rho == 0) {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("x", seq_len(p))
    return(X)
  }

  # For rho > 0, simulate one common latent factor and b independent noise
  # vectors for each block.
  G <- as.integer(p / b)
  X <- matrix(NA_real_, nrow = n, ncol = p)
  sr <- sqrt(rho)
  si <- sqrt(1 - rho)

  for (g in seq_len(G)) {
    cols <- ((g - 1L) * b + 1L):(g * b)
    z_common <- rnorm(n)
    z_noise <- matrix(rnorm(n * b), nrow = n, ncol = b)
    X[, cols] <- sr * z_common + si * z_noise
  }
  colnames(X) <- paste0("x", seq_len(p))
  return(X)
}

# PURPOSE
# Evaluate beta' Sigma beta without constructing the full p x p covariance
# matrix.
#
# WHY IT EXISTS
# This quadratic form is the population variance of X beta. Using the closed
# block formula is much faster and less memory-intensive than building Sigma,
# especially when p = 5,000.
#
# INPUTS
# beta: Coefficient vector.
# b: Block size.
# rho: Equicorrelation within every block.
#
# RETURNS
# The scalar population signal variance beta' Sigma beta.
quadform_block_equicorr <- function(beta, b, rho) {
  beta <- as.numeric(beta)
  b <- as.integer(b)
  G <- length(beta) / b
  V <- 0
  # Sum the analytical quadratic form separately over independent blocks.
  for (g in seq_len(G)) {
    idx <- ((g - 1L) * b + 1L):(g * b)
    x <- beta[idx]
    V <- V + (1 - rho) * sum(x^2) + rho * sum(x)^2
  }
  return(as.numeric(V))
}

# PURPOSE
# Construct the true coefficient vector for one setting and scale it to the
# requested population SNR.
#
# WHY IT EXISTS
# Signal support, magnitude, sign, and scaling jointly define the simulation.
# One implementation prevents the pilots from using subtly different truths
# and makes beta exactly reproducible from the setting and seed.
#
# STATISTICAL IDEA
# The first B_active blocks contain signal, and the first k variables within
# each of those blocks are active. Raw magnitudes are drawn independently from
# U(0.5, 1.5), but all active variables in the same block share one random
# sign. The entire raw vector is then multiplied by one constant so that
# beta' Sigma beta / sigma2 equals the setting-specific target SNR exactly.
#
# INPUTS
# setting_row: One row returned by make_setting() or make_all_settings().
# seed: Integer seed used only for signs and raw magnitudes.
# design: Shared design object.
#
# RETURNS
# A list containing the raw and scaled coefficient vectors, scaling details,
# deterministic support indices, block signs, coefficient-range diagnostics,
# and achieved theoretical SNR.
make_beta <- function(setting_row, seed, design = make_design()) {
  # A task-specific seed makes the coefficient draw exactly reproducible.
  set.seed(as.integer(seed))

  p <- as.integer(setting_row$p)
  b <- as.integer(setting_row$b)
  B_active <- as.integer(setting_row$B_active)
  k <- as.integer(setting_row$k)
  rho <- as.numeric(setting_row$rho)
  snr <- as.numeric(setting_row$snr)
  q_active <- as.integer(setting_row$q_active)

  beta0 <- numeric(p)
  active_blocks <- seq_len(B_active)
  active_variables <- integer(q_active)

  # Draw all block signs and raw magnitudes before filling the vector. One
  # common sign per block creates aligned effects within a block, while signs
  # may differ across active blocks.
  block_signs <- sample(c(-1L, 1L), size = B_active, replace = TRUE)
  magnitudes <- runif(
    n = q_active,
    min = design$beta_magnitude_min,
    max = design$beta_magnitude_max
  )

  # Fill the first k positions of every active block. pos tracks which segment
  # of the magnitude vector belongs to the current block.
  pos <- 1L
  for (ii in seq_len(B_active)) {
    g <- active_blocks[ii]
    block_cols <- ((g - 1L) * b + 1L):(g * b)
    active_in_block <- block_cols[seq_len(k)]
    take <- pos:(pos + k - 1L)

    active_variables[take] <- active_in_block
    beta0[active_in_block] <- block_signs[ii] * magnitudes[take]
    pos <- pos + k
  }

  # V0 is the population signal variance before SNR scaling.
  V0 <- quadform_block_equicorr(beta0, b, rho)
  if (!is.finite(V0) || V0 <= 0) {
    stop("make_beta(): the unscaled signal variance is not positive and finite.")
  }

  # Multiplying beta0 by this factor gives
  # Var(X beta) / sigma2 = snr.
  beta_scale <- sqrt(snr * design$sigma2 / V0)
  beta <- beta_scale * beta0

  if (any(!is.finite(beta[active_variables]))) {
    stop("make_beta(): a scaled active coefficient is not finite.")
  }

  achieved_snr <- quadform_block_equicorr(beta, b, rho) / design$sigma2

  return(list(
    beta0 = beta0,
    beta = beta,
    V0 = V0,
    beta_scale = beta_scale,
    active_blocks = as.integer(active_blocks),
    active_variables = as.integer(active_variables),
    block_signs = as.integer(block_signs),
    min_abs_beta0_active = min(abs(beta0[active_variables])),
    max_abs_beta0_active = max(abs(beta0[active_variables])),
    min_abs_beta_active = min(abs(beta[active_variables])),
    max_abs_beta_active = max(abs(beta[active_variables])),
    achieved_snr_theoretical = achieved_snr
  ))
}

# PURPOSE
# Generate one independent training sample and one independent test sample from
# the Gaussian linear data-generating mechanism for a fixed coefficient vector.
#
# WHY IT EXISTS
# Every method in a repetition must use exactly the same data, while different
# repetitions must be independently reproducible. Generating both samples once
# here ensures that model comparisons are paired fairly.
#
# INPUTS
# setting_row: One simulation setting, including n, p, rho, b, and n_test.
# beta_object: Output from make_beta(), beta_object$beta is the true vector.
# seed: Integer seed for predictors and Gaussian noise in both samples.
# design: Shared design object supplying sigma2.
#
# RETURNS
# A list with X_train, y_train, X_test, y_test and detailed generation times
# used by the runtime pilot.
simulate_data <- function(setting_row, beta_object, seed,
                          design = make_design()) {
  set.seed(as.integer(seed))

  t_total <- proc.time()[["elapsed"]]

  # Generate training predictors, noise, and response separately so their
  # contributions to runtime remain visible.
  tx <- proc.time()[["elapsed"]]
  X_train <- simulate_X_blockcor(
    setting_row$n, setting_row$p, setting_row$b, setting_row$rho
  )
  time_X_train_sec <- proc.time()[["elapsed"]] - tx

  te <- proc.time()[["elapsed"]]
  eps_train <- rnorm(setting_row$n, sd = sqrt(design$sigma2))
  time_eps_train_sec <- proc.time()[["elapsed"]] - te

  ty <- proc.time()[["elapsed"]]
  y_train <- as.numeric(X_train %*% beta_object$beta + eps_train)
  time_y_train_sec <- proc.time()[["elapsed"]] - ty

  # The test sample is independent of the training sample but uses the same
  # population covariance, coefficient vector, and noise variance.
  tx2 <- proc.time()[["elapsed"]]
  X_test <- simulate_X_blockcor(
    setting_row$n_test, setting_row$p, setting_row$b, setting_row$rho
  )
  time_X_test_sec <- proc.time()[["elapsed"]] - tx2

  te2 <- proc.time()[["elapsed"]]
  eps_test <- rnorm(setting_row$n_test, sd = sqrt(design$sigma2))
  time_eps_test_sec <- proc.time()[["elapsed"]] - te2

  ty2 <- proc.time()[["elapsed"]]
  y_test <- as.numeric(X_test %*% beta_object$beta + eps_test)
  time_y_test_sec <- proc.time()[["elapsed"]] - ty2

  return(list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test,
    time_X_train_sec = time_X_train_sec,
    time_eps_train_sec = time_eps_train_sec,
    time_y_train_sec = time_y_train_sec,
    time_X_test_sec = time_X_test_sec,
    time_eps_test_sec = time_eps_test_sec,
    time_y_test_sec = time_y_test_sec,
    time_data_sec = proc.time()[["elapsed"]] - t_total
  ))
}

###############################################################
# D. CROSS-VALIDATION, MODELS, AND PERFORMANCE MEASURES
###############################################################

# PURPOSE
# Create a balanced, reproducible fold assignment for cv.glmnet().
#
# WHY IT EXISTS
# Ridge, LASSO, and both Elastic Net variants must use identical folds within
# a repetition. Otherwise fold randomness would contaminate their paired
# performance comparison.
#
# INPUTS
# n_train: Number of training observations.
# nfolds: Number of cross-validation folds.
# seed: Integer seed controlling only the fold permutation.
#
# RETURNS
# Integer vector of length n_train with fold labels 1, ..., nfolds.
make_cv_foldid <- function(n_train,
                           nfolds = make_design()$glmnet_nfolds,
                           seed) {
  n_train <- as.integer(n_train)
  nfolds <- as.integer(nfolds)
  set.seed(as.integer(seed))
  return(as.integer(sample(
    rep(seq_len(nfolds), length.out = n_train),
    size = n_train,
    replace = FALSE
  )))
}

# PURPOSE
# Compute prediction mean squared error (MSE) and mean absolute error (MAE).
#
# WHY IT EXISTS
# Both prediction measures must be calculated from the same independent test
# observations and the same prediction-error vector. Keeping the calculation
# in one helper prevents model-specific code from applying different formulas.
#
# INPUTS
# y: Observed test outcomes.
# prediction: Model predictions in the same order.
#
# RETURNS
# A list with scalar mse and mae values. MSE is the primary performance metric,
# MAE is retained as a secondary robustness measure.
compute_prediction_metrics <- function(y, prediction) {
  if (length(y) != length(prediction)) {
    stop("The observed and predicted outcome vectors have different lengths.")
  }
  # Work with the common prediction-error vector so both metrics are based on
  # exactly the same test observations.
  err <- as.numeric(y) - as.numeric(prediction)
  mse <- mean(err^2)
  mae <- mean(abs(err))
  if (!is.finite(mse) || mse < 0 || !is.finite(mae) || mae < 0) {
    stop("A prediction metric is not finite or is negative.")
  }
  return(list(mse = mse, mae = mae))
}

# PURPOSE
# Measure squared coefficient-estimation error overall and separately on the
# true active and inactive coordinates.
#
# WHY IT EXISTS
# Prediction is primary, but coefficient error helps explain whether a method's
# predictive behavior is driven by shrinkage on signals or noise variables.
#
# INPUTS
# beta_hat: Estimated p-vector without the fitted intercept.
# beta_true: True p-vector used to generate the data.
# active_variables: Integer indices of nonzero true coefficients.
#
# RETURNS
# A list with overall, active-set, and inactive-set coefficient MSE.
compute_coefficient_mse <- function(beta_hat, beta_true, active_variables) {
  beta_hat <- as.numeric(beta_hat)
  beta_true <- as.numeric(beta_true)
  if (length(beta_hat) != length(beta_true)) {
    stop("The estimated and true coefficient vectors have different lengths.")
  }
  p <- length(beta_true)
  active <- sort(unique(as.integer(active_variables)))
  inactive <- setdiff(seq_len(p), active)
  # Compute coordinate-wise squared errors once, then average the relevant
  # subsets without refitting or transforming the model.
  sq <- (beta_hat - beta_true)^2

  return(list(
    coefficient_mse_overall = mean(sq),
    coefficient_mse_active = if (length(active) > 0) mean(sq[active]) else NA_real_,
    coefficient_mse_inactive = if (length(inactive) > 0) mean(sq[inactive]) else NA_real_
  ))
}

# PURPOSE
# Fit one penalized Gaussian regression model with cross-validation, predict on
# the independent test set, and calculate prediction and coefficient metrics.
#
# WHY IT EXISTS
# Ridge, LASSO, and both Elastic Net variants differ only in alpha. One shared
# wrapper guarantees identical glmnet options, folds, prediction rules, and
# timing definitions across methods.
#
# INPUTS
# dat: Output from simulate_data().
# alpha: glmnet mixing parameter (0 Ridge, 1 LASSO, intermediate Elastic Net).
# foldid: Shared fold vector returned by make_cv_foldid().
# beta_true: Optional true coefficient vector for coefficient MSE.
# active_variables: Optional active-coordinate indices.
# design: Shared model settings.
#
# RETURNS
# A list containing prediction mean squared error, prediction mean absolute
# error, fitted coefficients, the two penalty choices, coefficient mean squared
# errors, and separate timing components.
fit_glmnet <- function(dat, alpha, foldid, beta_true = NULL,
                       active_variables = NULL,
                       design = make_design()) {
  # Fit and select lambda using the supplied common folds. glmnet standardizes
  # predictors internally and fits an unpenalized intercept, as specified in
  # the Registered Report. The population intercept in the data-generating
  # mechanism is zero. Nevertheless, the fitted model retains an intercept
  # because a finite training sample need not be centered exactly.
  # parallel = FALSE is used because batchtools manages parallelism at the
  # task level.
  fit_time <- system.time({
    cv_fit <- glmnet::cv.glmnet(
      x = dat$X_train,
      y = dat$y_train,
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
  })

  # Evaluate the selected lambda.min model on the independent test matrix.
  pred_time <- system.time({
    prediction <- as.numeric(stats::predict(
      cv_fit,
      newx = dat$X_test,
      s = design$glmnet_prediction_s,
      type = "response"
    ))
  })

  # Time prediction and coefficient metrics separately for the runtime pilot.
  metric_time <- system.time({
    pm <- compute_prediction_metrics(dat$y_test, prediction)
  })

  coef_time <- system.time({
    cc <- as.matrix(stats::coef(cv_fit, s = design$glmnet_prediction_s))
    intercept_hat <- as.numeric(cc[1L, 1L])
    beta_hat <- as.numeric(cc[-1L, 1L])

    cm <- if (!is.null(beta_true) && !is.null(active_variables)) {
      compute_coefficient_mse(beta_hat, beta_true, active_variables)
    } else {
      list(
        coefficient_mse_overall = NA_real_,
        coefficient_mse_active = NA_real_,
        coefficient_mse_inactive = NA_real_
      )
    }
  })

  lambda_min <- as.numeric(cv_fit$lambda.min)
  lambda_1se <- as.numeric(cv_fit$lambda.1se)
  return(c(
    list(
      mse = pm$mse,
      mae = pm$mae,
      intercept_hat = intercept_hat,
      beta_hat = beta_hat,
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      fit_time_sec = as.numeric(fit_time["elapsed"]),
      pred_time_sec = as.numeric(pred_time["elapsed"]),
      prediction_metric_time_sec = as.numeric(metric_time["elapsed"]),
      coefficient_metric_time_sec = as.numeric(coef_time["elapsed"])
    ),
    cm
  ))
}

# PURPOSE
# Fit the Random Forest comparator used only by the blinded precision pilot.
#
# WHY IT EXISTS
# The precision pilot needs a paired variability pattern without fitting Ridge
# before the confirmatory simulation. Random Forest supplies that pattern. Its
# mean performance is not interpreted as Ridge performance.
#
# INPUTS
# dat: Output from simulate_data().
# p: Number of predictors, used for the fixed floor(sqrt(p)) mtry choice in
#   this pilot.
# seed: Integer seed for the forest.
# design: Tree count and minimum-node-size settings.
#
# RETURNS
# A list with prediction mean squared error, prediction mean absolute error,
# the number of candidate predictors considered at each split, and fitting,
# prediction, and metric times.
fit_rf <- function(dat, p, seed, design = make_design()) {
  p <- as.integer(p)
  mtry <- max(1L, as.integer(floor(sqrt(p))))

  # One ranger thread avoids nested parallelism inside batchtools workers.
  fit_time <- system.time({
    rf_fit <- ranger::ranger(
      x = dat$X_train,
      y = dat$y_train,
      num.trees = as.integer(design$rf_num_trees),
      mtry = mtry,
      min.node.size = as.integer(design$rf_min_node_size),
      seed = as.integer(seed),
      num.threads = 1L,
      write.forest = TRUE
    )
  })

  pred_time <- system.time({
    prediction <- as.numeric(stats::predict(
      rf_fit,
      data = dat$X_test,
      num.threads = 1L
    )$predictions)
  })

  metric_time <- system.time({
    pm <- compute_prediction_metrics(dat$y_test, prediction)
  })

  return(list(
    mse = pm$mse,
    mae = pm$mae,
    mtry = mtry,
    fit_time_sec = as.numeric(fit_time["elapsed"]),
    pred_time_sec = as.numeric(pred_time["elapsed"]),
    prediction_metric_time_sec = as.numeric(metric_time["elapsed"])
  ))
}

###############################################################
# E. PRIMARY RATIO-OF-MEANS ESTIMAND
###############################################################

# PURPOSE
# Estimate the relative performance contrast
#   theta = 1 - mean(comparator MSE) / mean(baseline MSE)
# and its paired Monte Carlo uncertainty.
#
# WHY IT EXISTS
# The study estimand is a ratio of method-specific mean prediction mean squared
# errors, not the mean of repetition-specific ratios. Because both methods are
# evaluated on the same repetitions, the influence-value calculation retains
# their covariance and gives the appropriate Monte Carlo standard error. This
# first-order calculation is the delta method for a ratio metric described by
# Deng et al. (2018).
#
# INPUTS
# baseline_mse: Repetition-level prediction mean squared error of the
#   denominator method, which is LASSO here.
# comparator_mse: Paired repetition-level prediction mean squared error of the
#   comparison method.
# conf_level: Two-sided normal-approximation confidence level.
#
# RETURNS
# A one-row tibble with sample size, method means, theta, influence-value
# standard deviation tau, Monte Carlo standard error, confidence interval, and
# paired correlation.
ratio_of_means_summary <- function(baseline_mse,
                                   comparator_mse,
                                   conf_level = 0.95) {
  # Retain only complete finite pairs, pairing must never be broken.
  ok <- is.finite(baseline_mse) & is.finite(comparator_mse)
  L <- as.numeric(baseline_mse[ok])
  C <- as.numeric(comparator_mse[ok])
  n <- length(L)
  if (n < 2L) {
    return(tibble(
      n = n,
      mean_baseline = if (n > 0) mean(L) else NA_real_,
      mean_comparator = if (n > 0) mean(C) else NA_real_,
      theta = NA_real_,
      tau_influence = NA_real_,
      mcse = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      correlation = NA_real_
    ))
  }

  mu_L <- mean(L)
  mu_C <- mean(C)
  if (!is.finite(mu_L) || mu_L <= 0) {
    stop("The mean baseline MSE must be positive and finite.")
  }

  # Calculate the ratio only after each method has been averaged.
  theta <- 1 - mu_C / mu_L

  # Paired influence value for g(mu_L, mu_C) = 1 - mu_C / mu_L.
  psi <-
    (mu_C / mu_L^2) * (L - mu_L) -
    (1 / mu_L) * (C - mu_C)

  # tau is the between-repetition standard deviation of the influence values.
  # Dividing it by the square root of n gives the Monte Carlo standard error of
  # the estimated theta.
  tau <- stats::sd(psi)
  mcse <- tau / sqrt(n)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  return(tibble(
    n = n,
    mean_baseline = mu_L,
    mean_comparator = mu_C,
    theta = theta,
    tau_influence = tau,
    mcse = mcse,
    ci_lower = theta - z * mcse,
    ci_upper = theta + z * mcse,
    correlation = stats::cor(L, C)
  ))
}

# PURPOSE
# Convert an estimated influence-value standard deviation into the required
# repetition count.
#
# WHY IT EXISTS
# If the Monte Carlo standard error is approximately tau divided by the square
# root of R, rearranging the target inequality gives R >= (tau/target)^2.
# Keeping the formula here makes the precision-planning rule explicit and
# applies the prespecified minimum in one place. Morris et al. (2019) and Siepe
# et al. (2024) motivate planning and reporting Monte Carlo uncertainty in
# simulation studies. The equation used here follows by directly rearranging
# the study's Monte Carlo standard-error approximation. Here tau does not denote
# an ordinary outcome standard deviation. It is the standard deviation of the
# paired influence value for the ratio-of-means estimand, so the repetition
# calculation targets the exact quantity reported by the study.
#
# INPUTS
# tau: Estimated standard deviation of the paired influence values.
# mcse_target: Desired maximum Monte Carlo standard error.
# R_min: Prespecified lower bound on repetitions.
#
# RETURNS
# Integer max(R_min, ceiling((tau / mcse_target)^2)), or NA for invalid tau.
required_R_from_tau <- function(tau,
                                mcse_target = make_design()$mcse_target,
                                R_min = make_design()$R_min) {
  if (!is.finite(tau) || tau < 0) return(NA_integer_)
  return(as.integer(max(R_min, ceiling((tau / mcse_target)^2))))
}

###############################################################
# F. TASKS AND RESULT TEMPLATES
###############################################################

# PURPOSE
# Expand the setting table into one row per setting and Monte Carlo repetition.
#
# WHY IT EXISTS
# Both serial execution and batchtools need an explicit task table. Storing
# setting_id, run_id, and seed in every row makes scheduling independent of the
# order in which workers finish.
#
# INPUTS
# settings: Output from make_all_settings().
# R: Number of repetitions per setting.
# seed_base: Integer offset defining an independent seed family for the pilot.
#
# RETURNS
# A tibble with all setting columns plus task_id, run_id, and a deterministic
# task seed. The seed formula reserves four component seeds per repetition and
# prevents their numeric reuse across repetitions and settings.
make_tasks <- function(settings, R, seed_base) {
  if (length(R) != 1L || !is.finite(R) || R != floor(R) ||
    R < 1 || R > 9999) {
    stop("R must be one integer from 1 to 9999.")
  }
  if (nrow(settings) < 1L || anyDuplicated(settings$setting_id)) {
    stop("settings must contain distinct setting IDs.")
  }
  R <- as.integer(R)
  # Validate the complete range before allocating or converting any seed.
  task_seeds <- make_task_seeds(
    rep(settings$setting_id, each = R),
    rep(seq_len(R), times = nrow(settings)), seed_base
  )
  rows <- vector("list", nrow(settings) * R)
  counter <- 1L

  # Copy each setting R times. task_id is a storage/scheduler identifier,
  # setting_id and run_id retain the statistical meaning.
  for (i in seq_len(nrow(settings))) {
    s <- settings[i, , drop = FALSE]
    for (r in seq_len(R)) {
      one <- s
      one$run_id <- r
      one$seed <- task_seeds[counter]
      one$task_id <- counter
      rows[[counter]] <- one
      counter <- counter + 1L
    }
  }
  return(bind_rows(rows) %>% relocate(task_id, setting_id, run_id, seed))
}

# PURPOSE
# Initialize the columns common to precision and runtime task results.
#
# WHY IT EXISTS
# A task may fail after only some steps have completed. Starting with explicit
# NA values guarantees that even failed tasks return a rectangular row with
# setting metadata, the failed step, and the error message.
#
# INPUTS
# task_row: One row from make_tasks().
#
# RETURNS
# A one-row tibble containing identifiers, design factors, placeholders for
# coefficient diagnostics and timings, and ok = FALSE until completion.
base_result_template <- function(task_row) {
  out <- tibble(
    task_id = as.integer(task_row$task_id),
    setting_id = as.integer(task_row$setting_id),
    run_id = as.integer(task_row$run_id),
    seed = as.integer(task_row$seed),
    n = as.integer(task_row$n),
    p = as.integer(task_row$p),
    p_over_n = as.numeric(task_row$p_over_n),
    rho = as.numeric(task_row$rho),
    a = as.numeric(task_row$a),
    k = as.integer(task_row$k),
    snr = as.numeric(task_row$snr),
    b = as.integer(task_row$b),
    G = as.integer(task_row$G),
    B_active = as.integer(task_row$B_active),
    q_active = as.integer(task_row$q_active),
    s_realized = as.numeric(task_row$s_realized),
    n_test = as.integer(task_row$n_test),
    V0 = NA_real_,
    beta_scale = NA_real_,
    min_abs_beta0_active = NA_real_,
    min_abs_beta_active = NA_real_,
    max_abs_beta_active = NA_real_,
    active_blocks = list(integer(0)),
    active_variables = list(integer(0)),
    block_signs = list(integer(0)),
    time_beta_sec = NA_real_,
    time_data_sec = NA_real_,
    time_total_sec = NA_real_,
    ok = FALSE,
    failed_step = NA_character_,
    error_message = NA_character_
  )
  return(out)
}

# PURPOSE
# Copy coefficient-generation metadata into an initialized result row.
#
# WHY IT EXISTS
# Both task types record the same compact data-generating diagnostics.
# Centralizing these assignments prevents the precision and runtime workers
# from storing different metadata. It also avoids saving the full,
# reproducible beta vector.
#
# INPUTS
# res: One-row tibble from base_result_template().
# beta_obj: Output from make_beta().
#
# RETURNS
# The same result row with SNR-scaling diagnostics and support/sign lists
# filled. The full p-vector is deliberately not stored because it is exactly
# reproducible from the setting and seed and would make result files much
# larger.
fill_beta_info <- function(res, beta_obj) {
  res$V0 <- beta_obj$V0
  res$beta_scale <- beta_obj$beta_scale
  res$min_abs_beta0_active <- beta_obj$min_abs_beta0_active
  res$min_abs_beta_active <- beta_obj$min_abs_beta_active
  res$max_abs_beta_active <- beta_obj$max_abs_beta_active

  res$active_blocks <- list(beta_obj$active_blocks)
  res$active_variables <- list(beta_obj$active_variables)
  res$block_signs <- list(beta_obj$block_signs)
  return(res)
}

# PURPOSE
# Execute one blinded precision-pilot repetition.
#
# WHY IT EXISTS
# The worker deliberately fits LASSO and Random Forest but not Ridge. Their
# paired prediction mean squared error variation is later transformed into
# conservative planning scenarios without examining the confirmatory
# Ridge-to-LASSO comparison.
#
# INPUTS
# task_row: One task produced by make_tasks().
# design: Shared design object.
#
# RETURNS
# One internal result row with the paired LASSO and Random Forest prediction
# mean squared errors needed for precision planning, model timings,
# data-generating metadata, and an explicit success or failure record. These
# individual prediction errors are not written to the public result files.
run_one_sd_task <- function(task_row, design = make_design()) {
  start <- proc.time()[["elapsed"]]
  task_row <- tibble::as_tibble(task_row)
  res <- base_result_template(task_row) %>%
    mutate(
      mse_lasso = NA_real_, mse_rf = NA_real_,
      lambda_lasso_min = NA_real_,
      time_lasso_fit_sec = NA_real_, time_lasso_pred_sec = NA_real_,
      time_rf_fit_sec = NA_real_, time_rf_pred_sec = NA_real_
    )
  step <- "initialization"

  # Each stochastic component receives a fixed offset from the task seed, so
  # changing one step does not silently change all later random draws.
  res <- tryCatch({
    # 1. Construct and SNR-scale the true coefficients.
    step <- "beta"
    tb <- proc.time()[["elapsed"]]
    beta_obj <- make_beta(task_row, task_row$seed + 1L, design)
    res$time_beta_sec <- proc.time()[["elapsed"]] - tb
    res <- fill_beta_info(res, beta_obj)

    # 2. Generate the paired training and independent test samples.
    step <- "data"
    dat <- simulate_data(task_row, beta_obj, task_row$seed + 2L, design)
    res$time_data_sec <- dat$time_data_sec

    # 3. Create folds once. Only LASSO uses them here, but this matches the
    # glmnet procedure used in the runtime pilot and the planned study.
    step <- "foldid"
    foldid <- make_cv_foldid(
      nrow(dat$X_train), design$glmnet_nfolds, task_row$seed + 3L
    )

    # 4. Fit LASSO and evaluate it on the common test set.
    step <- "lasso"
    lasso <- fit_glmnet(
      dat, alpha = 1, foldid = foldid,
      beta_true = beta_obj$beta,
      active_variables = beta_obj$active_variables,
      design = design
    )

    # 5. Fit the blinded Random Forest comparator on the same data.
    step <- "rf"
    rf <- fit_rf(dat, task_row$p, task_row$seed + 4L, design)

    res$mse_lasso <- lasso$mse
    res$mse_rf <- rf$mse
    res$lambda_lasso_min <- lasso$lambda_min
    res$time_lasso_fit_sec <- lasso$fit_time_sec
    res$time_lasso_pred_sec <- lasso$pred_time_sec
    res$time_rf_fit_sec <- rf$fit_time_sec
    res$time_rf_pred_sec <- rf$pred_time_sec
    # A task is successful only after all results have been copied.
    res$ok <- TRUE
    res
  }, error = function(e) {
    # Return the partial row instead of stopping the whole worker collection.
    # The summary function refuses to analyze any failed task.
    res$failed_step <- step
    res$error_message <- conditionMessage(e)
    return(res)
  })

  res$time_total_sec <- proc.time()[["elapsed"]] - start
  return(res)
}

# PURPOSE
# Execute and time one complete runtime-pilot repetition.
#
# WHY IT EXISTS
# Runtime planning must measure the same work later performed per repetition:
# coefficient/data generation, common folds, four cv.glmnet fits, predictions,
# and both prediction and coefficient metrics. Scientific method comparisons
# are not inferred from this function's output.
#
# INPUTS
# task_row: One task produced by make_tasks().
# design: Shared design object.
#
# RETURNS
# One result row with fine-grained timing components, selected lambdas,
# reproducibility metadata, and success/failure information.
run_one_runtime_task <- function(task_row, design = make_design()) {
  start <- proc.time()[["elapsed"]]
  task_row <- tibble::as_tibble(task_row)
  res <- base_result_template(task_row) %>%
    mutate(
      cv_fold_seed = as.integer(task_row$seed + 3L),
      cv_nfolds = as.integer(design$glmnet_nfolds),
      time_X_train_sec = NA_real_, time_eps_train_sec = NA_real_,
      time_y_train_sec = NA_real_, time_X_test_sec = NA_real_,
      time_eps_test_sec = NA_real_, time_y_test_sec = NA_real_,
      time_foldid_sec = NA_real_,
      time_ridge_fit_sec = NA_real_, time_ridge_pred_sec = NA_real_,
      time_ridge_prediction_metric_sec = NA_real_,
      time_ridge_coefficient_metric_sec = NA_real_,
      time_ridge_total_sec = NA_real_, time_ridge_wrapper_sec = NA_real_,
      time_lasso_fit_sec = NA_real_, time_lasso_pred_sec = NA_real_,
      time_lasso_prediction_metric_sec = NA_real_,
      time_lasso_coefficient_metric_sec = NA_real_,
      time_lasso_total_sec = NA_real_, time_lasso_wrapper_sec = NA_real_,
      time_en_fit_sec = NA_real_, time_en_pred_sec = NA_real_,
      time_en_prediction_metric_sec = NA_real_,
      time_en_coefficient_metric_sec = NA_real_,
      time_en_total_sec = NA_real_, time_en_wrapper_sec = NA_real_,
      time_en095_fit_sec = NA_real_, time_en095_pred_sec = NA_real_,
      time_en095_prediction_metric_sec = NA_real_,
      time_en095_coefficient_metric_sec = NA_real_,
      time_en095_total_sec = NA_real_, time_en095_wrapper_sec = NA_real_,
      time_prediction_metrics_sec = NA_real_,
      time_coefficient_metrics_sec = NA_real_,
      time_internal_unaccounted_sec = NA_real_,
      lambda_ridge_min = NA_real_, lambda_ridge_1se = NA_real_,
      lambda_lasso_min = NA_real_, lambda_lasso_1se = NA_real_,
      lambda_en_min = NA_real_, lambda_en_1se = NA_real_,
      lambda_en095_min = NA_real_, lambda_en095_1se = NA_real_
    )
  step <- "initialization"

  # As in the precision task, independent seed offsets isolate coefficient,
  # data, and fold randomness.
  res <- tryCatch({
    # 1. Generate the true coefficient vector and time that step.
    step <- "beta"
    tb <- proc.time()[["elapsed"]]
    beta_obj <- make_beta(task_row, task_row$seed + 1L, design)
    res$time_beta_sec <- proc.time()[["elapsed"]] - tb
    res <- fill_beta_info(res, beta_obj)

    # 2. Generate both data sets and retain the subcomponent timings returned
    # by simulate_data().
    step <- "data"
    dat <- simulate_data(task_row, beta_obj, task_row$seed + 2L, design)
    res$time_data_sec <- dat$time_data_sec
    res$time_X_train_sec <- dat$time_X_train_sec
    res$time_eps_train_sec <- dat$time_eps_train_sec
    res$time_y_train_sec <- dat$time_y_train_sec
    res$time_X_test_sec <- dat$time_X_test_sec
    res$time_eps_test_sec <- dat$time_eps_test_sec
    res$time_y_test_sec <- dat$time_y_test_sec

    # 3. Create the shared five-fold allocation once for all four methods.
    step <- "foldid"
    tf <- proc.time()[["elapsed"]]
    foldid <- make_cv_foldid(
      nrow(dat$X_train), design$glmnet_nfolds, task_row$seed + 3L
    )
    res$time_foldid_sec <- proc.time()[["elapsed"]] - tf

    # PURPOSE
    # Fit one alpha value and measure the complete wrapper time in addition to
    # the lower-level components returned by fit_glmnet().
    #
    # WHY IT EXISTS
    # The same timing definition must be used for all four methods. A local
    # function avoids duplicated stopwatch code while retaining access to this
    # task's data, folds, beta, active set, and design.
    #
    # INPUT
    # alpha: glmnet mixing value for the method being timed.
    #
    # RETURN VALUE
    # The list from fit_glmnet() with wrapper_elapsed_sec added.
    fit_one_timed <- function(alpha) {
      tm <- proc.time()[["elapsed"]]
      out <- fit_glmnet(
        dat, alpha, foldid,
        beta_true = beta_obj$beta,
        active_variables = beta_obj$active_variables,
        design = design
      )
      out$wrapper_elapsed_sec <- proc.time()[["elapsed"]] - tm
      return(out)
    }

    # 4. Fit all four methods on the same generated data and cross-validation
    # folds.
    step <- "ridge"
    ridge <- fit_one_timed(0)
    step <- "lasso"
    lasso <- fit_one_timed(1)
    step <- "elastic_net_1over3"
    en <- fit_one_timed(design$en_alpha)
    step <- "elastic_net_095"
    en095 <- fit_one_timed(design$en_near_lasso_alpha)

    # PURPOSE
    # Copy one model's named timing fields into the common task-result row.
    #
    # WHY IT EXISTS
    # Four methods expose identical timing components. Centralizing the mapping
    # prevents method-specific column mistakes and calculates the same wrapper
    # remainder for every fit.
    #
    # INPUTS
    # res: The task-result row to update.
    # prefix: Column-name prefix identifying ridge, lasso, en, or en095.
    # x: Output from fit_one_timed().
    #
    # RETURN VALUE
    # The updated task-result row. The wrapper field is total wrapper time
    # minus explicitly timed components.
    fill_model_times <- function(res, prefix, x) {
      res[[paste0("time_", prefix, "_fit_sec")]] <- x$fit_time_sec
      res[[paste0("time_", prefix, "_pred_sec")]] <- x$pred_time_sec
      res[[paste0("time_", prefix, "_prediction_metric_sec")]] <- x$prediction_metric_time_sec
      res[[paste0("time_", prefix, "_coefficient_metric_sec")]] <- x$coefficient_metric_time_sec
      component_sum <- x$fit_time_sec + x$pred_time_sec +
        x$prediction_metric_time_sec + x$coefficient_metric_time_sec
      res[[paste0("time_", prefix, "_total_sec")]] <- x$wrapper_elapsed_sec
      res[[paste0("time_", prefix, "_wrapper_sec")]] <-
        max(0, x$wrapper_elapsed_sec - component_sum)
      return(res)
    }

    res <- fill_model_times(res, "ridge", ridge)
    res <- fill_model_times(res, "lasso", lasso)
    res <- fill_model_times(res, "en", en)
    res <- fill_model_times(res, "en095", en095)

    # 5. Aggregate metric-computation time across the four methods.
    res$time_prediction_metrics_sec <- sum(c(
      ridge$prediction_metric_time_sec,
      lasso$prediction_metric_time_sec,
      en$prediction_metric_time_sec,
      en095$prediction_metric_time_sec
    ))
    res$time_coefficient_metrics_sec <- sum(c(
      ridge$coefficient_metric_time_sec,
      lasso$coefficient_metric_time_sec,
      en$coefficient_metric_time_sec,
      en095$coefficient_metric_time_sec
    ))

    # 6. Preserve the selected tuning parameters for reproducibility. The
    # runtime pilot does not interpret their relative scientific performance.
    res$lambda_ridge_min <- ridge$lambda_min
    res$lambda_ridge_1se <- ridge$lambda_1se
    res$lambda_lasso_min <- lasso$lambda_min
    res$lambda_lasso_1se <- lasso$lambda_1se
    res$lambda_en_min <- en$lambda_min
    res$lambda_en_1se <- en$lambda_1se
    res$lambda_en095_min <- en095$lambda_min
    res$lambda_en095_1se <- en095$lambda_1se

    res$ok <- TRUE
    res
  }, error = function(e) {
    # Keep a structured failure record so a worker error cannot disappear in
    # the scheduler logs. The summary refuses to use failed tasks.
    res$failed_step <- step
    res$error_message <- conditionMessage(e)
    return(res)
  })

  res$time_total_sec <- proc.time()[["elapsed"]] - start

  # Compare the complete task stopwatch with the sum of explicitly timed
  # components. The non-negative remainder includes R-level orchestration.
  accounted <- sum(c(
    res$time_beta_sec,
    res$time_data_sec,
    res$time_foldid_sec,
    res$time_ridge_total_sec,
    res$time_lasso_total_sec,
    res$time_en_total_sec,
    res$time_en095_total_sec
  ), na.rm = TRUE)
  res$time_internal_unaccounted_sec <- max(0, res$time_total_sec - accounted)
  return(res)
}

###############################################################
# G. DETAILED RUNTIME AND SCHEDULER SUMMARIES
###############################################################
# These helpers live in the core so the runtime tests can evaluate them
# without accidentally starting the full runtime pilot.
###############################################################

# PURPOSE
# Convert batchtools' job table into a task-linked scheduler timing table.
#
# WHY IT EXISTS
# Internal task stopwatches do not include queueing and scheduler overhead.
# When reliable batchtools timestamps are available, this table connects those
# external times to the corresponding statistical task_id.
#
# INPUTS
# task_run: The list returned by run_tasks(), including registry and results.
#
# RETURNS
# A tibble with job/task identifiers, submitted/started/done timestamps,
# queueing time, worker-running time, scheduler-reported memory when available,
# and estimated worker overhead. An empty tibble is returned for serial
# execution or an empty registry.
extract_job_runtime_detail <- function(task_run) {
  if (is.null(task_run$registry)) return(tibble())

  reg <- task_run$registry
  job_table <- batchtools::getJobTable(reg = reg)
  if (nrow(job_table) == 0L) return(tibble())

  # batchtools has used two names for the stored parameter column across
  # versions. Select whichever is present and recover task_id from it.
  pars_column <- if ("pars" %in% names(job_table)) {
    "pars"
  } else if ("job.pars" %in% names(job_table)) {
    "job.pars"
  } else {
    NA_character_
  }

  task_ids <- if (!is.na(pars_column)) {
    vapply(job_table[[pars_column]], function(x) {
      return(if (!is.null(x$task_row$task_id)) {
        as.integer(x$task_row$task_id)
      } else {
        NA_integer_
      })
    }, integer(1))
  } else {
    rep(NA_integer_, nrow(job_table))
  }

  # Copy timestamps once so the duration logic below is easy to read.
  submitted <- job_table$submitted
  started <- job_table$started
  done <- job_table$done

  # Some batchtools versions provide durations directly. Otherwise they are
  # reconstructed from the submitted, started, and completed timestamps.
  queued_sec <- if ("time.queued" %in% names(job_table)) {
    as.numeric(job_table$time.queued)
  } else {
    as.numeric(difftime(started, submitted, units = "secs"))
  }
  running_sec <- if ("time.running" %in% names(job_table)) {
    as.numeric(job_table$time.running)
  } else {
    as.numeric(difftime(done, started, units = "secs"))
  }

  return(tibble(
    job_id = as.integer(job_table$job.id),
    task_id = task_ids,
    submitted = submitted,
    started = started,
    done = done,
    time_queued_sec = queued_sec,
    time_running_sec = running_sec,
    time_submission_to_done_sec = queued_sec + running_sec,
    memory_used = if ("mem.used" %in% names(job_table)) {
      job_table$mem.used
    } else {
      NA
    },
    batchtools_error = if ("error" %in% names(job_table)) job_table$error else NA
  ) %>%
    left_join(
      task_run$run_results %>%
        select(task_id, setting_id, ok, time_total_sec),
      by = "task_id"
    ) %>%
    mutate(
      worker_overhead_sec = pmax(0, time_running_sec - time_total_sec)
    ))
}

# PURPOSE
# Summarize runtime-task measurements by setting and quantify scheduler
# overhead when complete timestamps are available.
#
# WHY IT EXISTS
# The later projection needs setting-specific mean task times because large-p
# settings have different runtimes from small-p settings. Scheduler overhead
# is reported separately
# and used only with complete timestamps and one recorded submission on a
# matching worker setup.
#
# INPUTS
# run_results: Task-level rows returned by run_one_runtime_task().
# settings: Full setting table used for the pilot.
# job_detail: Output from extract_job_runtime_detail().
# R_runtime: Repetitions per setting in the runtime pilot, stored as metadata.
# n_cores: Number of workers used for the pilot.
# design: Shared design object.
# execution_metadata: Persisted runtime execution history returned by run_tasks().
#
# RETURNS
# A list with by_setting (mean timing components plus design factors) and
# global (overall runtime, timestamp coverage, makespan, and overhead factor).
summarise_runtime_detail <- function(run_results, settings, job_detail,
                                     R_runtime, n_cores,
                                     design = make_design(),
                                     execution_metadata = NULL) {
  if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1 ||
      n_cores != as.integer(n_cores)) stop("n_cores must be one positive integer.")
  n_cores <- as.integer(n_cores)
  R_runtime <- as.integer(R_runtime)
  if (anyDuplicated(run_results$task_id) > 0L) {
    stop("The runtime results contain duplicate task_id values.")
  }
  if (any(!run_results$ok)) {
    stop("The runtime pilot contains failed tasks.")
  }

  # 1. Average each timing component within its setting. This preserves the
  # factorial structure for the later repetition-count projection.
  by_setting <- run_results %>%
    group_by(setting_id) %>%
    summarise(
      n_tasks = n(),
      n_success = sum(ok),
      n_failed = sum(!ok),
      mean_time_total_sec = safe_mean(time_total_sec),
      median_time_total_sec = safe_median(time_total_sec),
      mean_time_beta_sec = safe_mean(time_beta_sec),
      mean_time_data_sec = safe_mean(time_data_sec),
      mean_time_X_train_sec = safe_mean(time_X_train_sec),
      mean_time_X_test_sec = safe_mean(time_X_test_sec),
      mean_time_eps_train_sec = safe_mean(time_eps_train_sec),
      mean_time_eps_test_sec = safe_mean(time_eps_test_sec),
      mean_time_y_train_sec = safe_mean(time_y_train_sec),
      mean_time_y_test_sec = safe_mean(time_y_test_sec),
      mean_time_foldid_sec = safe_mean(time_foldid_sec),
      mean_time_ridge_fit_sec = safe_mean(time_ridge_fit_sec),
      mean_time_lasso_fit_sec = safe_mean(time_lasso_fit_sec),
      mean_time_en_fit_sec = safe_mean(time_en_fit_sec),
      mean_time_en095_fit_sec = safe_mean(time_en095_fit_sec),
      mean_time_ridge_pred_sec = safe_mean(time_ridge_pred_sec),
      mean_time_lasso_pred_sec = safe_mean(time_lasso_pred_sec),
      mean_time_en_pred_sec = safe_mean(time_en_pred_sec),
      mean_time_en095_pred_sec = safe_mean(time_en095_pred_sec),
      mean_time_prediction_metrics_sec = safe_mean(time_prediction_metrics_sec),
      mean_time_coefficient_metrics_sec = safe_mean(time_coefficient_metrics_sec),
      mean_time_internal_unaccounted_sec = safe_mean(time_internal_unaccounted_sec),
      .groups = "drop"
    ) %>%
    left_join(settings, by = "setting_id")

  # 2. Check the three scheduler timestamps for each job before
  # calculating scheduler overhead.
  timestamp_complete <- nrow(job_detail) > 0L &&
    nrow(job_detail) == nrow(run_results) &&
    identical(sort(as.integer(job_detail$task_id)),
              sort(as.integer(run_results$task_id))) &&
    all(!is.na(job_detail$submitted)) &&
    all(!is.na(job_detail$started)) &&
    all(!is.na(job_detail$done)) &&
    all(job_detail$submitted <= job_detail$started) &&
    all(job_detail$started <= job_detail$done)

  # A resumed benchmark mixes computing time with time between sessions.
  # Only one recorded submission on the matching worker setup supplies an
  # overhead factor. Internal task times remain available after interruption.
  metadata_complete <- is.list(execution_metadata) &&
    identical(execution_metadata$schema, "runtime_execution_v1") &&
    length(execution_metadata$submission_count) == 1L &&
    is.finite(execution_metadata$submission_count) &&
    identical(execution_metadata$profile$n_cores, n_cores) &&
    is.character(execution_metadata$profile$cluster_function) &&
    length(execution_metadata$profile$cluster_function) == 1L &&
    length(execution_metadata$profile$host) == 4L &&
    all(!is.na(execution_metadata$profile$host)) &&
    is.logical(execution_metadata$unknown_history) &&
    length(execution_metadata$unknown_history) == 1L &&
    !is.na(execution_metadata$unknown_history)
  continuous_benchmark <- metadata_complete &&
    execution_metadata$submission_count == 1L &&
    !execution_metadata$unknown_history
  overhead_valid <- timestamp_complete && continuous_benchmark
  overhead_note <- if (!timestamp_complete) {
    "Overhead unavailable: scheduler timestamps or task coverage are incomplete or invalid."
  } else if (!metadata_complete) {
    "Overhead unavailable: execution history or matching worker configuration is unverified."
  } else if (!continuous_benchmark) {
    "Overhead unavailable: benchmark was resumed or has unknown prior execution history."
  } else {
    "Overhead based on complete timestamps from one recorded submission on the matching worker setup."
  }

  submitted_coverage <- if (nrow(job_detail) > 0L) {
    mean(!is.na(job_detail$submitted))
  } else {
    NA_real_
  }

  # 3. With complete timestamps, compare observed wall-clock makespan with the
  # ideal task-time sum divided by the number of workers.
  job_global <- if (timestamp_complete) {
    start0 <- min(job_detail$submitted)
    end0 <- max(job_detail$done)
    makespan_sec <- as.numeric(difftime(end0, start0, units = "secs"))
    internal_sum <- sum(run_results$time_total_sec, na.rm = TRUE)
    ideal_sec <- internal_sum / n_cores

    # Scheduler overhead is used only with complete timestamps. Small timing
    # deviations below one are not interpreted as a speedup.
    raw_overhead <- if (overhead_valid && is.finite(makespan_sec) && ideal_sec > 0) {
      makespan_sec / ideal_sec
    } else {
      NA_real_
    }
    overhead_factor <- if (is.finite(raw_overhead)) max(1, raw_overhead) else NA_real_
    efficiency <- if (is.finite(overhead_factor)) 1 / overhead_factor else NA_real_

    tibble(
      runtime_timestamps_complete = TRUE,
      submitted_timestamp_coverage = submitted_coverage,
      pilot_makespan_sec = makespan_sec,
      pilot_makespan_hours = makespan_sec / 3600,
      ideal_parallel_sec = ideal_sec,
      observed_parallel_efficiency = efficiency,
      observed_overhead_factor = overhead_factor,
      mean_queue_sec = safe_mean(job_detail$time_queued_sec),
      median_queue_sec = safe_median(job_detail$time_queued_sec),
      mean_worker_running_sec = safe_mean(job_detail$time_running_sec),
      mean_worker_overhead_sec = safe_mean(job_detail$worker_overhead_sec),
      runtime_note = overhead_note
    )
  } else {
    # Keep internal timing summaries when scheduler timestamps are incomplete,
    # but do not manufacture an overhead adjustment from partial information.
    tibble(
      runtime_timestamps_complete = FALSE,
      submitted_timestamp_coverage = submitted_coverage,
      pilot_makespan_sec = NA_real_,
      pilot_makespan_hours = NA_real_,
      ideal_parallel_sec = sum(run_results$time_total_sec, na.rm = TRUE) / n_cores,
      observed_parallel_efficiency = NA_real_,
      observed_overhead_factor = NA_real_,
      mean_queue_sec = safe_mean(job_detail$time_queued_sec),
      median_queue_sec = safe_median(job_detail$time_queued_sec),
      mean_worker_running_sec = safe_mean(job_detail$time_running_sec),
      mean_worker_overhead_sec = safe_mean(job_detail$worker_overhead_sec),
      runtime_note = overhead_note
    )
  }

  # 4. Combine task-level summaries and scheduler-level summaries in one row.
  global <- tibble(
    design_name = design$design_name,
    R_runtime = as.integer(R_runtime),
    pilot_n_cores = n_cores,
    runtime_execution_metadata_complete = metadata_complete,
    runtime_submission_count = if (metadata_complete) {
      as.integer(execution_metadata$submission_count)
    } else NA_integer_,
    runtime_overhead_valid = overhead_valid,
    n_settings = nrow(settings),
    n_tasks = nrow(run_results),
    n_success = sum(run_results$ok),
    n_failed = sum(!run_results$ok),
    mean_time_total_sec = safe_mean(run_results$time_total_sec),
    median_time_total_sec = safe_median(run_results$time_total_sec),
    sum_internal_task_sec = sum(run_results$time_total_sec, na.rm = TRUE),
    n_test = design$n_test,
    created_at = Sys.time()
  ) %>%
    bind_cols(job_global)

  return(list(by_setting = by_setting, global = global))
}

# PURPOSE
# Project runtime-pilot measurements to a planned study with a specified
# number of repetitions per setting and n_cores parallel workers.
#
# WHY IT EXISTS
# Per-setting runtime differs substantially, so a projection should preserve
# those means rather than multiply one overall average. The function also
# separates an always-available ideal parallel estimate from an optional
# scheduler-adjusted estimate that is used only with complete timestamps.
# The calculation is
#   projected serial seconds
#     = repetitions per setting * sum(setting-specific mean task seconds),
#   ideal parallel seconds = projected serial seconds / number of workers,
#   adjusted seconds = ideal parallel seconds * observed overhead factor.
# The last line requires complete submission/start/completion timestamps,
# one recorded submission and the same planned worker count. Paused or
# unverified benchmarks provide internal-time projections only.
#
# INPUTS
# runtime_by_setting: Setting-level table from summarise_runtime_detail().
# runtime_global: Global table from summarise_runtime_detail().
# repetitions_per_setting: Planned repetitions for every design setting.
# n_cores: Planned number of parallel workers.
#
# RETURNS
# A one-row tibble with serial, ideal-parallel, and—only when justified by
# a verified continuous benchmark—scheduler-adjusted runtime projections.
project_planned_study_runtime <- function(runtime_by_setting, runtime_global,
                                          repetitions_per_setting,
                                          n_cores = 8L) {
  if (length(repetitions_per_setting) != 1L ||
      !is.finite(repetitions_per_setting) || repetitions_per_setting < 1 ||
      repetitions_per_setting != as.integer(repetitions_per_setting)) {
    stop("repetitions_per_setting must be one positive integer.")
  }
  if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1 ||
      n_cores != as.integer(n_cores)) stop("n_cores must be one positive integer.")
  repetitions_per_setting <- as.integer(repetitions_per_setting)
  n_cores <- as.integer(n_cores)

  if (any(!is.finite(runtime_by_setting$mean_time_total_sec))) {
    stop("At least one setting-level mean runtime is not finite.")
  }

  # Step 1: multiply every setting's observed mean task time by the requested
  # number of future repetitions, then sum across settings. This preserves the
  # large runtime differences between low-dimensional and high-dimensional
  # settings instead of replacing them with one unweighted overall mean.
  serial_sec <- sum(
    runtime_by_setting$mean_time_total_sec * repetitions_per_setting
  )
  ideal_sec <- serial_sec / n_cores

  timestamp_complete <- isTRUE(
    runtime_global[["runtime_timestamps_complete"]][1]
  )
  overhead_factor <- runtime_global[["observed_overhead_factor"]][1]

  # Step 2: decide whether the scheduler-based adjustment is defensible. The
  # ideal projection always remains available when timestamp metadata are
  # incomplete, but an adjusted value is not manufactured from partial data.
  matching_workers <- length(runtime_global[["pilot_n_cores"]]) == 1L &&
    isTRUE(runtime_global[["pilot_n_cores"]][1L] == n_cores)
  use_overhead <- timestamp_complete &&
    isTRUE(runtime_global[["runtime_overhead_valid"]][1L]) && matching_workers &&
    length(overhead_factor) == 1L &&
    is.finite(overhead_factor) &&
    overhead_factor >= 1

  adjusted_sec <- if (use_overhead) {
    ideal_sec * overhead_factor
  } else {
    NA_real_
  }

  # Step 3: return hours and days for interpretation while retaining the
  # assumptions needed to reconstruct the calculation.
  return(tibble(
    repetitions_per_setting = repetitions_per_setting,
    n_settings = nrow(runtime_by_setting),
    total_tasks = nrow(runtime_by_setting) * repetitions_per_setting,
    n_cores = n_cores,
    projected_serial_hours = serial_sec / 3600,
    projected_ideal_parallel_hours = ideal_sec / 3600,
    projected_ideal_parallel_days = ideal_sec / 86400,
    runtime_timestamps_complete = timestamp_complete,
    runtime_overhead_used = use_overhead,
    observed_overhead_factor_from_pilot =
      if (use_overhead) overhead_factor else NA_real_,
    projected_adjusted_hours = adjusted_sec / 3600,
    projected_adjusted_days = adjusted_sec / 86400,
    projection_note = if (use_overhead) {
      "Adjusted projection uses a verified continuous benchmark with the same worker count."
    } else if (!matching_workers && length(runtime_global[["pilot_n_cores"]]) == 1L) {
      "Ideal projection only: the requested worker count differs from the benchmark."
    } else {
      paste("Ideal projection only:",
            if (length(runtime_global[["runtime_note"]]) == 1L) runtime_global[["runtime_note"]]
            else "a continuous matching benchmark is not verified.")
    }
  ))
}

###############################################################
# H. BATCHTOOLS AND SERIAL EXECUTION
###############################################################

# PURPOSE
# Attach the requested execution backend to an existing batchtools registry.
#
# WHY IT EXISTS
# The same registry code supports a one-worker diagnostic run and the study's
# multi-worker socket run. Configuring the backend in one place prevents the
# task orchestration from depending on an implicit machine default.
#
# INPUTS
# reg: Writable batchtools registry.
# n_cores: Number of task workers.
# cluster_function: "socket" for the study server, or another explicitly
# supported batchtools backend.
#
# RETURNS
# The registry after saving the selected cluster function.
set_batchtools_cluster_function <- function(reg, n_cores,
                                            cluster_function = "socket") {
  n_cores <- as.integer(n_cores)

  # A requested socket backend always creates real socket workers, including
  # a one-job smoke test. An in-process substitute would miss worker-loading
  # failures. Other one-worker backends retain the original interactive path.
  if (cluster_function == "socket") {
    if (!requireNamespace("snow", quietly = TRUE)) {
      stop("The socket backend requires package 'snow'. Run run/00_install_packages.R.")
    }
    reg$cluster.functions <- batchtools::makeClusterFunctionsSocket(n_cores)
  } else if (n_cores == 1L) {
    reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive(external = FALSE)
  } else if (cluster_function == "multicore") {
    reg$cluster.functions <- batchtools::makeClusterFunctionsMulticore(n_cores)
  } else if (cluster_function == "interactive_external") {
    reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive(external = TRUE)
  } else {
    stop("Unknown cluster_function: ", cluster_function)
  }

  batchtools::saveRegistry(reg = reg)
  return(reg)
}


# PURPOSE
# Determine whether a directory already contains a batchtools registry.
#
# WHY IT EXISTS
# Resuming requires loading an existing registry, whereas an execution without
# a registry must create one. Testing for batchtools' registry.rds file
# distinguishes those paths.
#
# INPUTS
# registry_dir: Proposed registry directory.
#
# RETURNS
# TRUE when both the directory and registry.rds exist, otherwise FALSE.
registry_exists <- function(registry_dir) {
  return(dir.exists(registry_dir) && file.exists(file.path(registry_dir, "registry.rds")))
}


# PURPOSE
# Create a stable MD5 fingerprint of an arbitrary R object.
#
# WHY IT EXISTS
# Resume safety requires comparing the exact task keys and design with the
# registry's specification. Serializing to a temporary R data file gives
# a stable byte representation for that comparison.
#
# INPUTS
# x: Any serializable R object.
#
# RETURNS
# One MD5 character string. The temporary file is always deleted on exit.
object_md5 <- function(x) {
  tmp <- tempfile("pilot_obj_", fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  return(unname(tools::md5sum(tmp))[1])
}


# PURPOSE
# Fingerprint a worker function's formal arguments and body.
#
# WHY IT EXISTS
# A registry created with a different worker implementation must not be resumed,
# because its completed jobs and the requested jobs would use different code.
#
# INPUTS
# fun: R function used as the task worker.
#
# RETURNS
# One MD5 character string based on deparsed formals and function body.
function_md5 <- function(fun) {

  tmp <- tempfile("pilot_fun_")
  on.exit(unlink(tmp), add = TRUE)

  txt <- c(
    paste(deparse(formals(fun), width.cutoff = 500L), collapse = "\n"),
    paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
  )
  writeLines(txt, tmp, useBytes = TRUE)
  return(unname(tools::md5sum(tmp))[1])
}


# PURPOSE
# Identify which components differ between a saved and current run signature.
#
# WHY IT EXISTS
# A registry may be resumed only when its statistical design, task definition,
# seeds, and worker code agree with the requested execution. Naming every
# mismatching component produces an interpretable error instead of silently
# combining incompatible results.
#
# INPUTS
# saved: Signature read from an existing registry.
# current: Signature constructed for the requested run.
#
# RETURNS
# A character vector of mismatching field names. An empty vector means that the
# two signatures are compatible.
signature_mismatches <- function(saved, current) {
  nms <- union(names(saved), names(current))
  return(nms[!vapply(
    nms,
    function(nm) {
      return(nm %in% names(saved) &&
        nm %in% names(current) &&
        identical(saved[[nm]], current[[nm]]))
    },
    logical(1)
  )])
}


# PURPOSE
# Derive the lock-file path associated with one registry directory.
#
# WHY IT EXISTS
# The lock must be addressable before the registry itself exists. Storing it
# beside, rather than inside, the registry also protects creation and deletion.
#
# INPUTS
# registry_dir: Registry directory that requires exclusive write access.
#
# RETURNS
# Character path to the registry's companion lock file.
registry_lock_path <- function(registry_dir) {
  return(file.path(
    dirname(registry_dir),
    paste0(basename(registry_dir), ".pilot_registry_lock.rds")
  ))
}


# PURPOSE
# Return a stable name for the machine running the current R process.
#
# WHY IT EXISTS
# A registry lock can safely be declared stale only when its process identifier
# belongs to this same host. The host name is therefore stored with the PID.
#
# INPUTS
# None. The function queries Sys.info() and then the HOSTNAME environment value.
#
# RETURNS
# One character machine name, with "unknown" as the final fallback.
current_host <- function() {
  x <- Sys.info()[["nodename"]]
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    x <- Sys.getenv("HOSTNAME", unset = "unknown")
  }
  return(as.character(x))
}


# PURPOSE
# Determine whether a process recorded in a registry lock is still alive.
#
# WHY IT EXISTS
# A dead local writer leaves a stale lock that can be removed safely. A live or
# remote writer must not be disturbed because concurrent registry writes can
# corrupt or duplicate a long simulation run.
#
# INPUTS
# pid: Recorded process identifier.
# host: Recorded machine name.
#
# RETURNS
# TRUE when the local process can be confirmed alive, FALSE when confirmed
# dead, and NA when the status cannot safely be determined. A lock created on
# another host is one case in which the status is unknown.
pid_alive <- function(pid, host) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || pid < 1L) return(NA)
  if (!identical(as.character(host), current_host())) return(NA)

  if (.Platform$OS.type == "unix" && dir.exists("/proc")) {
    return(dir.exists(file.path("/proc", as.character(pid))))
  }

  if (identical(pid, Sys.getpid())) return(TRUE)
  return(NA)
}


# PURPOSE
# Acquire exclusive write access to a batchtools registry.
#
# WHY IT EXISTS
# Two R sessions writing the same registry can duplicate or corrupt work. The
# lock records host, process ID, and start time. A lock is removed automatically
# only when its local process can be confirmed dead.
#
# INPUTS
# registry_dir: Directory of the registry to protect.
#
# RETURNS
# Path to the newly created lock file. Stops if a live or unverifiable lock
# already exists.
acquire_registry_lock <- function(registry_dir) {
  dir.create(dirname(registry_dir), recursive = TRUE, showWarnings = FALSE)
  lock_file <- registry_lock_path(registry_dir)

  # Inspect an existing lock before writing a new one.
  if (file.exists(lock_file)) {
    old <- tryCatch(readRDS(lock_file), error = function(e) return(NULL))
    if (is.null(old)) {
      stop("The registry lock exists but cannot be read: ", lock_file)
    }

    alive <- pid_alive(old$pid, old$host)
    if (isFALSE(alive)) {
      message("Removing a stale registry lock from a terminated process.")
      unlink(lock_file)
    } else {
      stop(
        "The registry is already locked. Host=", old$host,
        ", PID=", old$pid, ", started=", as.character(old$started_at), ". ",
        "Do not use two writing R sessions at the same time."
      )
    }
  }

  # Persist enough information to diagnose the writer and distinguish a stale
  # local lock from a potentially live remote one.
  saveRDS(
    list(
      host = current_host(),
      pid = Sys.getpid(),
      started_at = Sys.time(),
      registry_dir = normalizePath(
        registry_dir, winslash = "/", mustWork = FALSE
      )
    ),
    lock_file
  )

  return(lock_file)
}


# PURPOSE
# Remove a registry lock after a safely completed or confirmed-idle run.
#
# WHY IT EXISTS
# acquire_registry_lock() blocks concurrent writers. Its companion cleanup is
# needed after successful completion and during on.exit() error handling.
#
# INPUTS
# lock_file: Exact path returned by acquire_registry_lock().
#
# RETURNS
# invisible(TRUE), allowing cleanup code to call it without printing output.
release_registry_lock <- function(lock_file) {
  if (length(lock_file) == 1L && file.exists(lock_file)) unlink(lock_file)
  return(invisible(TRUE))
}

# PURPOSE
# Execute a task table sequentially in the current R session.
#
# WHY IT EXISTS
# Serial execution is useful for smoke tests and debugging because errors are
# easy to reproduce without a scheduler. Full pilots use batchtools instead.
#
# INPUTS
# tasks: Task table from make_tasks().
# task_function: Worker such as run_one_sd_task().
# design: Shared design passed to every worker call.
# progress_every: Frequency of progress messages.
#
# RETURNS
# One tibble formed by row-binding all task-function results.
run_tasks_serial <- function(tasks, task_function,
                             design = make_design(), progress_every = 50L) {
  task_list <- split(tasks, seq_len(nrow(tasks)))
  out <- vector("list", length(task_list))

  for (i in seq_along(task_list)) {
    out[[i]] <- task_function(task_list[[i]], design)
    if (i %% progress_every == 0L || i == length(task_list)) {
      message("Completed: ", i, " / ", length(task_list))
    }
  }

  return(bind_rows(out))
}


# PURPOSE
# Describe every run component that must remain unchanged when a batchtools
# registry is resumed.
#
# WHY IT EXISTS
# Resuming after code, seeds, settings, or task keys change would mix results
# from different experiments. The signature guards only registry compatibility,
# it does not assert a fixed number of scientifically required runs.
#
# INPUTS
# tasks: Complete task table.
# task_function: Worker function mapped to jobs.
# design: Shared design object.
# source_files: R files sourced by batchtools workers.
# execution_profile: Optional runtime worker/backend/host configuration.
#
# RETURNS
# A named list containing design/task metadata plus checksums of task keys,
# worker code, and source files. This is operational metadata, not a scientific
# test of how many tasks the design should contain.
make_task_run_signature <- function(tasks, task_function, design,
                                    source_files, execution_profile = NULL) {
  # The shared helper is part of the executed code, even when workers receive
  # only the core path explicitly. Its checksum also protects a resumed run.
  source_files <- unique(c(source_files, reproducibility_script_file))
  return(list(
    signature_schema = "task_runner_signature",
    execution_profile = execution_profile,
    design = design,
    n_tasks = nrow(tasks),
    R = as.integer(max(tasks$run_id)),
    seed_min = as.integer(min(tasks$seed)),
    seed_max = as.integer(max(tasks$seed)),
    task_keys_md5 = object_md5(
      as.data.frame(
        tasks %>%
          select(task_id, setting_id, run_id, seed) %>%
          arrange(task_id)
      )
    ),
    task_function_md5 = function_md5(task_function),
    source_md5 = unname(tools::md5sum(source_files)),
    settings = as.data.frame(
      tasks %>%
        select(setting_id, n, p_over_n, p, rho, a, k, snr, b, n_test) %>%
        distinct() %>%
        arrange(setting_id)
    )
  ))
}


# PURPOSE
# Create or safely resume a batchtools registry, submit unfinished jobs, wait
# for completion, and collect task results.
#
# WHY IT EXISTS
# Precision and runtime pilots are long enough to require resumable execution.
# This wrapper prevents accidental reuse of a registry with different code or
# settings, blocks simultaneous writers, and never silently ignores job errors.
#
# INPUTS
# tasks: Complete table from make_tasks().
# task_function: Worker function executed for each row.
# design: Shared design object passed to workers.
# registry_dir: Persistent registry directory.
# n_cores: Number of workers.
# cluster_function: batchtools backend name.
# resume: Reuse a compatible existing registry.
# overwrite_registry: Deliberately replace an existing registry.
# resubmit_errors: Deliberately reset previously failed jobs.
# wait: Must remain TRUE so results are checked before collection.
# source_files: Code files sourced in every worker process.
# runtime_benchmark: Record submission history and protect the timing setup.
#
# RETURNS
# A list containing the registry object/path, batchtools status, and combined
# task-level results.
run_tasks_batchtools <- function(tasks,
                                 task_function,
                                 design = make_design(),
                                 registry_dir,
                                 n_cores = 1L,
                                 cluster_function = design$batchtools_cluster_function,
                                 resume = TRUE,
                                 overwrite_registry = FALSE,
                                 resubmit_errors = FALSE,
                                 wait = TRUE,
                                 source_files = current_script_file,
                                 runtime_benchmark = FALSE) {
  # batchtools copies source files into its registry. Include this dependency
  # so the copied core can source its sibling in every fresh worker process.
  source_files <- unique(c(source_files, reproducibility_script_file))
  if (missing(registry_dir) || length(registry_dir) != 1L) {
    stop("registry_dir is missing.")
  }
  if (nrow(tasks) < 1L) stop("tasks is empty.")
  if (anyDuplicated(tasks$task_id) > 0L) stop("task_id is not unique.")
  if (!isTRUE(wait)) {
    stop("The protected pilot runner requires wait = TRUE.")
  }
  if (length(source_files) < 1L || any(!file.exists(source_files))) {
    stop("All source_files must exist before batchtools is started.")
  }

  # 1. Protect the registry before inspecting, creating, or deleting it.
  dir.create(dirname(registry_dir), recursive = TRUE, showWarnings = FALSE)
  lock_file <- acquire_registry_lock(registry_dir)
  reg <- NULL

  # Release the lock automatically only when no job remains on the system.
  on.exit({
    # Keep the lock while batchtools still sees jobs on the system.
    on_system_now <- if (!is.null(reg)) {
      tryCatch(nrow(batchtools::findOnSystem(reg = reg)), error = function(e) return(NA_integer_))
    } else {
      0L
    }

    if (!is.na(on_system_now) && on_system_now == 0L) {
      release_registry_lock(lock_file)
    } else {
      message("The registry lock remains because the system status is not safely idle.")
    }
  }, add = TRUE)

  # 2. Build the signature before deciding whether an existing registry can be
  # resumed safely.
  signature_file <- file.path(registry_dir, "run_signature.rds")
  execution_profile <- if (isTRUE(runtime_benchmark)) {
    list(n_cores = as.integer(n_cores), cluster_function = cluster_function,
         host = unname(Sys.info()[c("nodename", "sysname", "release", "machine")]))
  } else NULL
  current_signature <- make_task_run_signature(
    tasks, task_function, design, source_files,
    execution_profile = execution_profile
  )

  exists_now <- registry_exists(registry_dir)

  if (exists_now && !resume && !overwrite_registry) {
    stop(
      "The registry already exists. Change resume or overwrite_registry."
    )
  }

  if (exists_now && resume && !overwrite_registry) {
    # Resume only after an exact signature comparison.
    reg <- batchtools::loadRegistry(
      file.dir = registry_dir,
      writeable = TRUE,
      make.default = FALSE
    )

    if (!file.exists(signature_file)) stop("The registry has no run_signature.rds file.")
    saved_signature <- readRDS(signature_file)
    bad <- signature_mismatches(saved_signature, current_signature)

    if (length(bad) > 0L) {
      stop(
        "The registry does not match the current pilot run. Different fields: ",
        paste(bad, collapse = ", "), "."
      )
    }
  } else {
    # A new registry maps each task row to one independent batchtools job.
    if (dir.exists(registry_dir) && overwrite_registry) {
      unlink(registry_dir, recursive = TRUE, force = TRUE)
    }

    reg <- batchtools::makeRegistry(
      file.dir = registry_dir,
      packages = required_packages,
      source = source_files,
      seed = 1L,
      make.default = FALSE
    )
    reg <- set_batchtools_cluster_function(reg, n_cores, cluster_function)

    task_list <- split(tasks, seq_len(nrow(tasks)))
    batchtools::batchMap(
      fun = task_function,
      task_row = task_list,
      more.args = list(design = design),
      reg = reg
    )

    saveRDS(current_signature, signature_file)
    batchtools::saveRegistry(reg = reg)
  }

  # Persist timing history before submission, including interrupted attempts.
  # A later call that only collects completed jobs does not add a submission.
  runtime_execution_file <- file.path(registry_dir, "runtime_execution.rds")
  runtime_execution <- NULL
  if (isTRUE(runtime_benchmark)) {
    runtime_execution <- if (exists_now && !overwrite_registry &&
                             file.exists(runtime_execution_file)) {
      readRDS(runtime_execution_file)
    } else {
      list(schema = "runtime_execution_v1", submission_count = 0L,
           unknown_history = isTRUE(exists_now && !overwrite_registry),
           profile = execution_profile)
    }
    if (!identical(runtime_execution$profile, execution_profile) ||
        !identical(runtime_execution$schema, "runtime_execution_v1") ||
        length(runtime_execution$submission_count) != 1L ||
        !is.finite(runtime_execution$submission_count) ||
        runtime_execution$submission_count < 0L) {
      stop("The runtime execution history does not match the benchmark configuration.")
    }
  }

  # 3. Refresh registry state and refuse to reset jobs that a scheduler still
  # reports as running.
  reg <- set_batchtools_cluster_function(reg, n_cores, cluster_function)
  batchtools::syncRegistry(reg = reg)

  on_system_before <- batchtools::findOnSystem(reg = reg)
  if (nrow(on_system_before) > 0L) {
    stop(
      "The registry still reports ", nrow(on_system_before),
      " jobs on the system. No automatic reset or duplicate start is allowed. ",
      "Inspect the status again later."
    )
  }

  # 4. Previously failed jobs require an explicit resubmission decision.
  errors_before <- batchtools::findErrors(reg = reg)
  if (nrow(errors_before) > 0L && !isTRUE(resubmit_errors)) {
    stop(
      "The registry contains ", nrow(errors_before),
      " batchtools errors. Inspect the logs first or deliberately set ",
      "resubmit_errors = TRUE."
    )
  }

  if (nrow(errors_before) > 0L && isTRUE(resubmit_errors)) {
    batchtools::resetJobs(ids = errors_before, reg = reg)
  }

  # 5. Submit only unfinished jobs, already completed results remain intact.
  ids_to_submit <- batchtools::findNotDone(reg = reg)

  if (nrow(ids_to_submit) > 0L) {
    # On resume, reset only unfinished jobs, completed jobs remain untouched.
    if (exists_now) {
      batchtools::resetJobs(ids = ids_to_submit, reg = reg)
      batchtools::syncRegistry(reg = reg)
      ids_to_submit <- batchtools::findNotDone(reg = reg)
    }

    message(
      "Completed jobs: ", nrow(batchtools::findDone(reg = reg)),
      ", jobs to submit: ", nrow(ids_to_submit), "."
    )

    if (isTRUE(runtime_benchmark)) {
      runtime_execution$submission_count <- runtime_execution$submission_count + 1L
      saveRDS(runtime_execution, runtime_execution_file)
    }
    batchtools::submitJobs(ids = ids_to_submit, reg = reg)
    batchtools::waitForJobs(
      ids = ids_to_submit,
      timeout = Inf,
      stop.on.error = FALSE,
      stop.on.expire = FALSE,
      reg = reg
    )
  } else {
    message("All jobs in the registry are already complete.")
  }

  # 6. After waiting, synchronize again and stop if any job failed or remains
  # incomplete. Partial pilot outputs are never summarized as complete.
  batchtools::syncRegistry(reg = reg)

  errors_after <- batchtools::findErrors(reg = reg)
  not_done_after <- batchtools::findNotDone(reg = reg)

  if (nrow(errors_after) > 0L) {
    stop("The run finished with ", nrow(errors_after), " batchtools errors.")
  }
  if (nrow(not_done_after) > 0L) {
    stop("After the run, ", nrow(not_done_after), " jobs are still incomplete.")
  }

  # 7. Collect one result per completed job and combine the returned rows.
  status <- batchtools::getStatus(reg = reg)
  result_list <- batchtools::reduceResultsList(missing.val = NULL, reg = reg)
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  run_results <- if (length(result_list) > 0L) bind_rows(result_list) else tibble()

  release_registry_lock(lock_file)

  return(list(
    registry = reg,
    registry_dir = registry_dir,
    status = status,
    runtime_execution = runtime_execution,
    run_results = run_results
  ))
}


# PURPOSE
# Select serial or batchtools execution through one common interface.
#
# WHY IT EXISTS
# The same pilot logic is used for small local smoke runs and resumable parallel
# production runs. A common return structure keeps higher-level pilot functions
# independent of the selected execution mechanism.
#
# INPUTS
# The task table, worker function, design, and backend-specific options passed
# through to run_tasks_serial() or run_tasks_batchtools().
#
# RETURNS
# In both modes, a list with registry metadata, status, and run_results. The
# common shape keeps pilot launchers independent of the execution backend.
run_tasks <- function(tasks, task_function,
                      design = make_design(),
                      backend = c("batchtools", "serial"),
                      registry_dir = NULL,
                      n_cores = 1L,
                      cluster_function = design$batchtools_cluster_function,
                      resume = TRUE,
                      overwrite_registry = FALSE,
                      resubmit_errors = FALSE,
                      wait = TRUE,
                      source_files = current_script_file,
                      runtime_benchmark = FALSE) {
  backend <- match.arg(backend)

  if (backend == "serial") {
    return(list(
      registry = NULL,
      registry_dir = NA_character_,
      status = tibble(),
      run_results = run_tasks_serial(tasks, task_function, design)
    ))
  }

  return(run_tasks_batchtools(
    tasks = tasks,
    task_function = task_function,
    design = design,
    registry_dir = registry_dir,
    n_cores = n_cores,
    cluster_function = cluster_function,
    resume = resume,
    overwrite_registry = overwrite_registry,
    resubmit_errors = resubmit_errors,
    wait = wait,
    source_files = source_files,
    runtime_benchmark = runtime_benchmark
  ))
}

###############################################################
# I. USAGE
###############################################################
# Precision and runtime launchers are kept in separate scripts. This preserves
# one authoritative data-generating mechanism while keeping the pilots'
# scientific purposes distinct.
###############################################################
