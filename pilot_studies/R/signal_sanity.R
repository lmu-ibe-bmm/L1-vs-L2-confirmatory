# Check signal calibration and RF learning over all 192 Gaussian settings.

signal_sanity_version <- "1.2.0"

signal_sanity_required_core_functions <- c(
  "make_design", "make_all_settings", "make_beta", "simulate_data",
  "set_batchtools_cluster_function"
)

# Check that the shared DGM and socket functions are loaded.
# Inputs: none. Returns invisible TRUE or stops on a missing function.
signal_sanity_assert_core <- function() {
  missing <- signal_sanity_required_core_functions[
    !vapply(signal_sanity_required_core_functions, exists, logical(1), mode = "function", inherits = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "Required project function(s) are missing after sourcing pilot_core.R: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  return(invisible(TRUE))
}

# Validate one task field before converting its type.
# Inputs: x and its field name. Returns the single non-missing value.
signal_sanity_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x)) {
    stop("Expected one non-missing value for '", name, "'.", call. = FALSE)
  }
  return(x[[1L]])
}

# Use the shared Gaussian DGM for the signal check.
# Inputs: n_test validates the requested analysis size. Returns the native design.
# The single-task function selects those rows from the generated test sample.
signal_sanity_prepare_design <- function(n_test = 5000L) {
  signal_sanity_assert_core()
  n_test <- as.integer(n_test)
  if (length(n_test) != 1L || is.na(n_test) || n_test < 100L) {
    stop("n_test must be one integer >= 100.", call. = FALSE)
  }

  design <- make_design()
  if (!is.list(design)) {
    stop("make_design() did not return a list-like design object.", call. = FALSE)
  }

  if (is.null(design$sigma2) || length(design$sigma2) != 1L ||
    !is.finite(design$sigma2) || design$sigma2 <= 0) {
    stop("design$sigma2 must exist and be a single positive finite value.", call. = FALSE)
  }

  if (!is.null(design$rf_num_trees) && as.integer(design$rf_num_trees) != 200L) {
    warning(
      "design$rf_num_trees is ", design$rf_num_trees,
      ", but this diagnostic intentionally uses 200 trees to match the precision pilot."
    )
  }

  return(design)
}

# Build and validate all 192 settings with their existing IDs.
# Inputs: design and expected_n_settings. Returns the ordered full grid.
signal_sanity_make_settings <- function(design, expected_n_settings = 192L) {
  settings <- make_all_settings(design)
  settings <- as.data.frame(settings, stringsAsFactors = FALSE)

  if (nrow(settings) != as.integer(expected_n_settings)) {
    stop(
      "Expected ", expected_n_settings, " settings (the current full precision-pilot grid), but make_all_settings() returned ",
      nrow(settings), ". The study design appears to have changed, stopping rather than silently running a different grid.",
      call. = FALSE
    )
  }

  required <- c("n", "p_over_n", "p", "rho", "a", "k", "snr", "b")
  missing <- setdiff(required, names(settings))
  if (length(missing) > 0L) {
    stop("Settings are missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!"setting_id" %in% names(settings)) settings$setting_id <- seq_len(nrow(settings))
  settings$setting_id <- as.integer(settings$setting_id)
  if (anyNA(settings$setting_id) || anyDuplicated(settings$setting_id)) {
    stop("setting_id must be complete and unique.", call. = FALSE)
  }
  if (!identical(sort(settings$setting_id), seq_len(nrow(settings)))) {
    stop("Expected setting_id to be exactly 1,...,n_settings.", call. = FALSE)
  }

  expected_p <- as.integer(round(settings$n * settings$p_over_n))
  if (any(expected_p != as.integer(settings$p))) {
    bad <- which(expected_p != as.integer(settings$p))[1L]
    stop("Inconsistent n, p/n and p in setting_id ", settings$setting_id[bad], ".", call. = FALSE)
  }
  if (any(as.integer(settings$b) <= 0L) || any(as.integer(settings$p) %% as.integer(settings$b) != 0L)) {
    stop("All settings must have positive block size b and p divisible by b.", call. = FALSE)
  }

  return(settings[order(settings$setting_id), , drop = FALSE])
}

# Assign corrected seeds to every setting/repetition pair.
# Inputs: settings, R and seed_base. Returns task IDs, run IDs and base seeds.
signal_sanity_make_tasks <- function(settings, R = 50L, seed_base = 1300000000L) {

  if (length(R) != 1L || !is.finite(R) || R != floor(R) ||
    R < 1 || R > 9999) {
    stop("R must be one integer from 1 to 9999.", call. = FALSE)
  }
  if (nrow(settings) < 1L || anyDuplicated(settings$setting_id)) {
    stop("Settings must have distinct setting IDs.", call. = FALSE)
  }
  make_task_seeds(settings$setting_id, rep(R, nrow(settings)), seed_base)
  R <- as.integer(R)

  grid <- expand.grid(
    setting_id = as.integer(settings$setting_id),
    run_id = seq_len(R),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$setting_id, grid$run_id), , drop = FALSE]
  rownames(grid) <- NULL
  grid$task_id <- seq_len(nrow(grid))

  grid$seed <- make_task_seeds(grid$setting_id, grid$run_id, seed_base)

  components <- unlist(lapply(1:4, function(off) return(as.double(grid$seed) + off)), use.names = FALSE)
  if (anyDuplicated(components)) {
    stop("Internal error: component-seed collision detected.", call. = FALSE)
  }

  return(grid[, c("task_id", "setting_id", "run_id", "seed")])
}

# Check the population signal variance under block equicorrelation.
# Inputs: beta, rho and block size b. Returns beta transpose Sigma beta.
signal_sanity_analytic_signal_variance <- function(beta, rho, b) {
  beta <- as.numeric(beta)
  rho <- as.numeric(rho)
  b <- as.integer(b)
  if (length(beta) < 1L || any(!is.finite(beta))) stop("beta must be finite and non-empty.", call. = FALSE)
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) stop("rho must be in [0,1).", call. = FALSE)
  if (length(b) != 1L || is.na(b) || b < 1L || length(beta) %% b != 0L) {
    stop("b must be a positive divisor of length(beta).", call. = FALSE)
  }

  block_id <- rep(seq_len(length(beta) / b), each = b)
  within_sq <- sum(beta^2)
  block_sums <- vapply(split(beta, block_id), sum, numeric(1))
  return((1 - rho) * within_sq + rho * sum(block_sums^2))
}

# Calculate test MSE from aligned outcomes y and predictions pred.
# Returns one mean squared prediction error.
signal_sanity_mse <- function(y, pred) return(mean((as.numeric(y) - as.numeric(pred))^2))

# Calculate test MAE from aligned outcomes y and predictions pred.
# Returns one mean absolute prediction error.
signal_sanity_mae <- function(y, pred) return(mean(abs(as.numeric(y) - as.numeric(pred))))

# Check population/sample SNR and whether RF improves on the training-mean baseline.
# Inputs: task_row, settings, design, rf_num_trees and analysis_n_test.
# Returns one row with SNR, MSE/MAE, relative gains, seeds and timings.
# The RF, baseline and known-coefficient predictions use the same test rows.
# Population SNR is analytic. Sample SNR uses the realized signal/noise variances.
# Failures return ok = FALSE with task identifiers and an error message.
signal_sanity_one_task <- function(task_row, settings, design, rf_num_trees = 200L, analysis_n_test = 5000L) {
  started <- proc.time()[[3L]]

  task_id <- setting_id <- run_id <- seed <- NA_integer_
  analysis_n_test <- as.integer(analysis_n_test)
  if (length(analysis_n_test) != 1L || is.na(analysis_n_test) || analysis_n_test < 100L) {
    stop("analysis_n_test must be one integer >= 100.", call. = FALSE)
  }
  return(tryCatch({

    task_row <- as.data.frame(task_row, stringsAsFactors = FALSE)
    if (nrow(task_row) != 1L) stop("task_row must contain exactly one row.")

    task_id <- as.integer(signal_sanity_scalar(task_row$task_id, "task_id"))
    setting_id <- as.integer(signal_sanity_scalar(task_row$setting_id, "setting_id"))
    run_id <- as.integer(signal_sanity_scalar(task_row$run_id, "run_id"))
    seed <- as.integer(signal_sanity_scalar(task_row$seed, "seed"))

    idx <- match(setting_id, settings$setting_id)
    if (is.na(idx)) stop("Unknown setting_id: ", setting_id)
    setting <- settings[idx, , drop = FALSE]

    worker_row <- setting
    worker_row$task_id <- task_id
    worker_row$run_id <- run_id
    worker_row$seed <- seed

    beta_obj <- make_beta(setting_row = worker_row, seed = seed + 1L)
    if (is.null(beta_obj$beta)) stop("make_beta() result has no $beta component.")
    beta <- as.numeric(beta_obj$beta)

    dat <- simulate_data(
      setting_row = worker_row,
      beta_object = beta_obj,
      seed = seed + 2L,
      design = design
    )
    required_dat <- c("X_train", "y_train", "X_test", "y_test")
    missing_dat <- setdiff(required_dat, names(dat))
    if (length(missing_dat) > 0L) stop("simulate_data() result is missing: ", paste(missing_dat, collapse = ", "))

    X_train <- dat$X_train
    X_test_native <- dat$X_test
    y_train <- as.numeric(dat$y_train)
    y_test_native <- as.numeric(dat$y_test)
    p <- ncol(X_train)
    if (is.null(p) || p != length(beta) || ncol(X_test_native) != p) stop("X/beta dimension mismatch.")
    if (length(y_train) != nrow(X_train) || length(y_test_native) != nrow(X_test_native)) stop("X/y dimension mismatch.")

    native_n_test <- nrow(X_test_native)
    if (native_n_test < analysis_n_test) {
      stop(
        "Project-native test sample has only ", native_n_test,
        " rows, but this diagnostic requests ", analysis_n_test,
        ". Lower analysis_n_test or inspect the current project test-size design."
      )
    }

    analysis_idx <- seq_len(analysis_n_test)
    X_test <- X_test_native[analysis_idx, , drop = FALSE]
    y_test <- y_test_native[analysis_idx]

    rho <- as.numeric(setting$rho[[1L]])
    b <- as.integer(setting$b[[1L]])
    target_snr <- as.numeric(setting$snr[[1L]])
    sigma2 <- as.numeric(design$sigma2)
    if (length(target_snr) != 1L || !is.finite(target_snr) || target_snr <= 0) stop("Target SNR must be positive and finite.")

    signal_var_analytic <- signal_sanity_analytic_signal_variance(beta, rho = rho, b = b)
    analytic_snr <- signal_var_analytic / sigma2
    analytic_snr_rel_error <- analytic_snr / target_snr - 1

    eta_train <- drop(X_train %*% beta)
    eta_test <- drop(X_test %*% beta)
    eps_train <- y_train - eta_train
    eps_test <- y_test - eta_test

    signal_var_train <- stats::var(eta_train)
    signal_var_test <- stats::var(eta_test)
    noise_var_train <- stats::var(eps_train)
    noise_var_test <- stats::var(eps_test)
    empirical_snr_train <- signal_var_train / noise_var_train
    empirical_snr_test <- signal_var_test / noise_var_test

    pred_intercept <- rep(mean(y_train), length(y_test))
    pred_oracle <- eta_test

    rf_num_trees <- as.integer(rf_num_trees)
    if (rf_num_trees < 1L) stop("rf_num_trees must be positive.")
    rf_mtry <- max(1L, as.integer(floor(sqrt(p))))
    rf_min_node_size <- if (!is.null(design$rf_min_node_size)) as.integer(design$rf_min_node_size) else 5L
    if (length(rf_min_node_size) != 1L || is.na(rf_min_node_size) || rf_min_node_size < 1L) {
      stop("Invalid RF min.node.size in design.")
    }

    rf_seed <- seed + 4L
    rf_fit_start <- proc.time()[[3L]]
    rf_fit <- ranger::ranger(
      x = X_train,
      y = y_train,
      num.trees = rf_num_trees,
      mtry = rf_mtry,
      min.node.size = rf_min_node_size,
      num.threads = 1L,
      seed = rf_seed,
      write.forest = TRUE
    )
    rf_fit_sec <- proc.time()[[3L]] - rf_fit_start

    rf_pred_start <- proc.time()[[3L]]
    pred_rf <- predict(rf_fit, data = X_test, num.threads = 1L)$predictions
    rf_pred_sec <- proc.time()[[3L]] - rf_pred_start

    mse_intercept <- signal_sanity_mse(y_test, pred_intercept)
    mse_oracle <- signal_sanity_mse(y_test, pred_oracle)
    mse_rf <- signal_sanity_mse(y_test, pred_rf)
    mae_intercept <- signal_sanity_mae(y_test, pred_intercept)
    mae_oracle <- signal_sanity_mae(y_test, pred_oracle)
    mae_rf <- signal_sanity_mae(y_test, pred_rf)

    denom_mse <- mse_intercept - mse_oracle
    rf_gain <- 1 - mse_rf / mse_intercept
    oracle_gain <- 1 - mse_oracle / mse_intercept
    fraction_oracle_gain <- if (is.finite(denom_mse) && denom_mse > 0) {
      (mse_intercept - mse_rf) / denom_mse
    } else NA_real_
    rf_excess_fraction <- if (is.finite(denom_mse) && denom_mse > 0) {
      (mse_rf - mse_oracle) / denom_mse
    } else NA_real_

    denom_mae <- mae_intercept - mae_oracle
    rf_mae_gain <- 1 - mae_rf / mae_intercept
    oracle_mae_gain <- 1 - mae_oracle / mae_intercept
    fraction_oracle_mae_gain <- if (is.finite(denom_mae) && denom_mae > 0) {
      (mae_intercept - mae_rf) / denom_mae
    } else NA_real_

    result <- list(
      ok = TRUE,
      error_message = NA_character_,
      task_id = task_id,
      setting_id = setting_id,
      run_id = run_id,
      seed = seed,
      beta_seed = seed + 1L,
      data_seed = seed + 2L,
      rf_seed = rf_seed,
      n = as.integer(setting$n[[1L]]),
      p_over_n = as.numeric(setting$p_over_n[[1L]]),
      p = as.integer(setting$p[[1L]]),
      rho = rho,
      a = as.numeric(setting$a[[1L]]),
      k = as.integer(setting$k[[1L]]),
      snr = target_snr,
      b = b,
      n_test = as.integer(nrow(X_test)),
      n_test_generated = as.integer(native_n_test),
      sigma2 = sigma2,
      signal_var_analytic = signal_var_analytic,
      analytic_snr = analytic_snr,
      analytic_snr_rel_error = analytic_snr_rel_error,
      signal_var_train = signal_var_train,
      signal_var_test = signal_var_test,
      noise_var_train = noise_var_train,
      noise_var_test = noise_var_test,
      empirical_snr_train = empirical_snr_train,
      empirical_snr_test = empirical_snr_test,
      population_oracle_gain_reference = target_snr / (1 + target_snr),
      mse_intercept = mse_intercept,
      mse_oracle = mse_oracle,
      mse_rf = mse_rf,
      rf_mse_difference_vs_intercept = mse_intercept - mse_rf,
      oracle_mse_difference_vs_intercept = mse_intercept - mse_oracle,
      rf_gain_vs_intercept = rf_gain,
      oracle_gain_vs_intercept = oracle_gain,
      rf_fraction_oracle_gain = fraction_oracle_gain,
      rf_excess_risk_fraction = rf_excess_fraction,
      rf_beats_intercept = isTRUE(mse_rf < mse_intercept),
      oracle_beats_intercept = isTRUE(mse_oracle < mse_intercept),
      oracle_mse_over_sigma2 = mse_oracle / sigma2,
      intercept_mse_over_population_reference = mse_intercept / (sigma2 * (1 + target_snr)),
      mae_intercept = mae_intercept,
      mae_oracle = mae_oracle,
      mae_rf = mae_rf,
      rf_mae_gain_vs_intercept = rf_mae_gain,
      oracle_mae_gain_vs_intercept = oracle_mae_gain,
      rf_fraction_oracle_mae_gain = fraction_oracle_mae_gain,
      rf_num_trees = rf_num_trees,
      rf_mtry = rf_mtry,
      rf_min_node_size = rf_min_node_size,
      rf_fit_sec = rf_fit_sec,
      rf_pred_sec = rf_pred_sec,
      elapsed_sec = proc.time()[[3L]] - started
    )
    as.data.frame(result, stringsAsFactors = FALSE)

  }, error = function(e) {
    return(as.data.frame(list(
      ok = FALSE,
      error_message = conditionMessage(e),
      task_id = task_id,
      setting_id = setting_id,
      run_id = run_id,
      seed = seed,
      elapsed_sec = proc.time()[[3L]] - started
    ), stringsAsFactors = FALSE))
  }))
}

# Batchtools wrapper with the same inputs and output as signal_sanity_one_task().
# Delegates the calculation unchanged to the single-task function.
signal_sanity_bt_worker <- function(task_row, settings, design, rf_num_trees = 200L, analysis_n_test = 5000L) {
  return(signal_sanity_one_task(
    task_row = task_row,
    settings = settings,
    design = design,
    rf_num_trees = rf_num_trees,
    analysis_n_test = analysis_n_test
  ))
}

# Summarize finite repetition-level values with a descriptive t interval.
# Inputs: x and confidence level. Returns mean, Monte Carlo SE and bounds.
# The SE is sd(x) / sqrt(n). Fewer than two values give missing bounds.
signal_sanity_mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(c(mean = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  m <- mean(x)
  if (n == 1L) return(c(mean = m, se = NA_real_, lower = NA_real_, upper = NA_real_))
  se <- stats::sd(x) / sqrt(n)
  crit <- stats::qt((1 + level) / 2, df = n - 1L)
  return(c(mean = m, se = se, lower = m - crit * se, upper = m + crit * se))
}

# Aggregate successful repetitions within each setting.
# Inputs: results. Returns SNR/MSE checks, mean gains and descriptive t intervals.
# Positive RF gain means lower test MSE than the training-mean baseline.
signal_sanity_setting_summary <- function(results) {
  ok <- results[!is.na(results$ok) & results$ok, , drop = FALSE]
  if (nrow(ok) == 0L) stop("No successful task results to summarise.", call. = FALSE)

  split_res <- split(ok, ok$setting_id)

  rows <- lapply(split_res, function(d) {
    gain_ci <- signal_sanity_mean_ci(d$rf_gain_vs_intercept)
    fraction_ci <- signal_sanity_mean_ci(d$rf_fraction_oracle_gain)
    diff_ci <- signal_sanity_mean_ci(d$rf_mse_difference_vs_intercept)
    snr_ci <- signal_sanity_mean_ci(d$empirical_snr_test)

    return(data.frame(
      setting_id = as.integer(d$setting_id[1L]),
      n_success = nrow(d),
      n = as.integer(d$n[1L]),
      p_over_n = as.numeric(d$p_over_n[1L]),
      p = as.integer(d$p[1L]),
      rho = as.numeric(d$rho[1L]),
      a = as.numeric(d$a[1L]),
      k = as.integer(d$k[1L]),
      snr = as.numeric(d$snr[1L]),
      b = as.integer(d$b[1L]),
      n_test = as.integer(d$n_test[1L]),
      min_n_test_generated = min(as.integer(d$n_test_generated)),
      max_n_test_generated = max(as.integer(d$n_test_generated)),
      mean_analytic_snr = mean(d$analytic_snr),
      max_abs_analytic_snr_rel_error = max(abs(d$analytic_snr_rel_error)),
      mean_empirical_snr_test = snr_ci[["mean"]],
      empirical_snr_test_ci_lower = snr_ci[["lower"]],
      empirical_snr_test_ci_upper = snr_ci[["upper"]],
      mean_mse_intercept = mean(d$mse_intercept),
      sd_mse_intercept = stats::sd(d$mse_intercept),
      mean_mse_oracle = mean(d$mse_oracle),
      sd_mse_oracle = stats::sd(d$mse_oracle),
      mean_mse_rf = mean(d$mse_rf),
      sd_mse_rf = stats::sd(d$mse_rf),
      mean_rf_mse_difference_vs_intercept = diff_ci[["mean"]],
      rf_mse_difference_ci_lower = diff_ci[["lower"]],
      rf_mse_difference_ci_upper = diff_ci[["upper"]],
      mean_rf_gain_vs_intercept = gain_ci[["mean"]],
      median_rf_gain_vs_intercept = stats::median(d$rf_gain_vs_intercept),
      sd_rf_gain_vs_intercept = stats::sd(d$rf_gain_vs_intercept),
      rf_gain_ci_lower = gain_ci[["lower"]],
      rf_gain_ci_upper = gain_ci[["upper"]],
      q10_rf_gain_vs_intercept = as.numeric(stats::quantile(d$rf_gain_vs_intercept, 0.10, names = FALSE)),
      q90_rf_gain_vs_intercept = as.numeric(stats::quantile(d$rf_gain_vs_intercept, 0.90, names = FALSE)),
      prop_rf_beats_intercept = mean(d$rf_beats_intercept),
      mean_oracle_gain_vs_intercept = mean(d$oracle_gain_vs_intercept),
      mean_population_oracle_gain_reference = mean(d$population_oracle_gain_reference),
      mean_rf_fraction_oracle_gain = fraction_ci[["mean"]],
      median_rf_fraction_oracle_gain = stats::median(d$rf_fraction_oracle_gain, na.rm = TRUE),
      rf_fraction_oracle_gain_ci_lower = fraction_ci[["lower"]],
      rf_fraction_oracle_gain_ci_upper = fraction_ci[["upper"]],
      mean_rf_excess_risk_fraction = mean(d$rf_excess_risk_fraction, na.rm = TRUE),
      mean_oracle_mse_over_sigma2 = mean(d$oracle_mse_over_sigma2),
      mean_intercept_mse_over_population_reference = mean(d$intercept_mse_over_population_reference),
      mean_mae_intercept = mean(d$mae_intercept),
      mean_mae_oracle = mean(d$mae_oracle),
      mean_mae_rf = mean(d$mae_rf),
      mean_rf_mae_gain_vs_intercept = mean(d$rf_mae_gain_vs_intercept),
      mean_elapsed_sec = mean(d$elapsed_sec),
      stringsAsFactors = FALSE
    ))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out[order(out$setting_id), , drop = FALSE])
}

# Describe the distribution of setting-level RF gains within each SNR level.
# Inputs: setting_summary. Returns gain ranges and counts of positive intervals.
signal_sanity_by_snr <- function(setting_summary) {
  groups <- split(setting_summary, setting_summary$snr)

  rows <- lapply(groups, function(d) {
    return(data.frame(
      snr = as.numeric(d$snr[1L]),
      n_settings = nrow(d),
      min_mean_rf_gain_vs_intercept = min(d$mean_rf_gain_vs_intercept),
      q25_mean_rf_gain_vs_intercept = as.numeric(stats::quantile(d$mean_rf_gain_vs_intercept, 0.25, names = FALSE)),
      median_mean_rf_gain_vs_intercept = stats::median(d$mean_rf_gain_vs_intercept),
      q75_mean_rf_gain_vs_intercept = as.numeric(stats::quantile(d$mean_rf_gain_vs_intercept, 0.75, names = FALSE)),
      max_mean_rf_gain_vs_intercept = max(d$mean_rf_gain_vs_intercept),
      min_prop_rf_beats_intercept = min(d$prop_rf_beats_intercept),
      n_settings_rf_gain_ci_lower_gt_0 = sum(d$rf_gain_ci_lower > 0, na.rm = TRUE),
      n_settings_mean_rf_gain_le_0 = sum(d$mean_rf_gain_vs_intercept <= 0, na.rm = TRUE),
      min_mean_rf_fraction_oracle_gain = min(d$mean_rf_fraction_oracle_gain, na.rm = TRUE),
      median_mean_rf_fraction_oracle_gain = stats::median(d$mean_rf_fraction_oracle_gain, na.rm = TRUE),
      max_abs_analytic_snr_rel_error = max(d$max_abs_analytic_snr_rel_error),
      stringsAsFactors = FALSE
    ))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out[order(out$snr), , drop = FALSE])
}

# Summarize completeness, SNR accuracy and the least improved setting.
# Inputs: results, setting_summary and expected_jobs. Returns one overview row.
signal_sanity_overview <- function(results, setting_summary, expected_jobs) {
  n_ok <- sum(!is.na(results$ok) & results$ok)
  n_failed <- sum(is.na(results$ok) | !results$ok)
  worst_idx <- which.min(setting_summary$mean_rf_gain_vs_intercept)
  worst <- setting_summary[worst_idx, , drop = FALSE]

  return(data.frame(
    version = signal_sanity_version,
    expected_jobs = as.integer(expected_jobs),
    n_results = nrow(results),
    n_success = n_ok,
    n_failed = n_failed,
    n_settings = nrow(setting_summary),
    repetitions_per_setting = min(setting_summary$n_success),
    max_abs_analytic_snr_rel_error = max(setting_summary$max_abs_analytic_snr_rel_error),
    n_settings_mean_rf_gain_le_0 = sum(setting_summary$mean_rf_gain_vs_intercept <= 0, na.rm = TRUE),
    n_settings_rf_gain_ci_lower_gt_0 = sum(setting_summary$rf_gain_ci_lower > 0, na.rm = TRUE),
    min_mean_rf_gain_vs_intercept = min(setting_summary$mean_rf_gain_vs_intercept),
    median_mean_rf_gain_vs_intercept = stats::median(setting_summary$mean_rf_gain_vs_intercept),
    max_mean_rf_gain_vs_intercept = max(setting_summary$mean_rf_gain_vs_intercept),
    worst_setting_id = as.integer(worst$setting_id),
    worst_setting_n = as.integer(worst$n),
    worst_setting_p = as.integer(worst$p),
    worst_setting_rho = as.numeric(worst$rho),
    worst_setting_a = as.numeric(worst$a),
    worst_setting_k = as.integer(worst$k),
    worst_setting_snr = as.numeric(worst$snr),
    worst_setting_mean_rf_gain = as.numeric(worst$mean_rf_gain_vs_intercept),
    worst_setting_rf_gain_ci_lower = as.numeric(worst$rf_gain_ci_lower),
    worst_setting_rf_gain_ci_upper = as.numeric(worst$rf_gain_ci_upper),
    min_mean_rf_fraction_oracle_gain = min(setting_summary$mean_rf_fraction_oracle_gain, na.rm = TRUE),
    median_mean_rf_fraction_oracle_gain = stats::median(setting_summary$mean_rf_fraction_oracle_gain, na.rm = TRUE),
    min_prop_rf_beats_intercept = min(setting_summary$prop_rf_beats_intercept, na.rm = TRUE),
    analysis_note = "Diagnostic screen, effect magnitudes and CIs are descriptive, not new confirmatory hypotheses.",
    stringsAsFactors = FALSE
  ))
}

# Normalize batchtools job IDs returned as a table or vector.
# Inputs: x. Returns an integer vector, possibly empty.
signal_sanity_job_ids <- function(x) {
  if (is.null(x) || length(x) == 0L) return(integer())
  if (is.data.frame(x)) {
    if ("job.id" %in% names(x)) return(as.integer(x$job.id))
    if ("job.id" %in% colnames(x)) return(as.integer(x[, "job.id"]))
  }
  return(as.integer(x))
}

# Create a worker script that loads the shared core and signal module.
# Inputs: destination, project_root and source paths. Returns the written path.
signal_sanity_write_worker_bootstrap <- function(
    destination,
    project_root = ".",
    pilot_core_file = "R/pilot_core.R",
    signal_sanity_file = "R/signal_sanity.R"
) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  pilot_abs <- normalizePath(file.path(project_root, pilot_core_file), winslash = "/", mustWork = TRUE)
  sanity_abs <- normalizePath(file.path(project_root, signal_sanity_file), winslash = "/", mustWork = TRUE)

  destination <- path.expand(destination)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  q <- function(x) return(encodeString(x, quote = '"'))
  lines <- c(
    "# Auto-generated signal-sanity batchtools worker bootstrap.",
    "# It deliberately sources the real project files by absolute path.",
    'Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",',
    '           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")',
    paste0("source(", q(pilot_abs), ")"),
    paste0("source(", q(sanity_abs), ")")
  )
  writeLines(lines, con = destination, useBytes = TRUE)

  if (!file.exists(destination)) {
    stop("Failed to create batchtools worker bootstrap: ", destination, call. = FALSE)
  }
  return(normalizePath(destination, winslash = "/", mustWork = TRUE))
}

# Record the diagnostic options, settings and source hashes before reusing jobs.
# Inputs: run options, settings and source_files. Returns a registry signature.
signal_sanity_signature <- function(R, n_test, rf_num_trees, n_cores, seed_base, settings, source_files) {

  source_files <- unique(c(source_files, reproducibility_script_file))
  source_files <- normalizePath(source_files, winslash = "/", mustWork = TRUE)
  md5 <- as.character(tools::md5sum(source_files))
  names(md5) <- basename(source_files)
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(settings, tf, version = 2)
  settings_md5 <- unname(tools::md5sum(tf))
  return(list(
    version = signal_sanity_version,
    R = as.integer(R),
    n_test = as.integer(n_test),
    rf_num_trees = as.integer(rf_num_trees),
    seed_base = as.integer(seed_base),
    n_settings = nrow(settings),
    setting_id = as.integer(settings$setting_id),
    settings_md5 = settings_md5,
    source_md5 = md5,
    cluster_function = "socket",
    n_cores = as.integer(n_cores)
  ))
}

# Combine worker rows while retaining short error records.
# Inputs: result_list. Returns one data frame with absent fields filled by NA.
signal_sanity_bind_results <- function(result_list) {
  good <- Filter(Negate(is.null), result_list)
  if (length(good) == 0L) return(data.frame())

  rows <- lapply(good, function(x) {
    if (is.data.frame(x)) return(x)
    return(as.data.frame(x, stringsAsFactors = FALSE))
  })
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))

  rows <- lapply(rows, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    return(d[, all_names, drop = FALSE])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}

# Run one real signal task in a temporary socket registry.
# Inputs: settings, design, task_row and reduced test/RF options.
# Returns the checked result and removes the temporary registry.
signal_sanity_batchtools_smoke_test <- function(
    settings,
    design,
    task_row,
    analysis_n_test = 500L,
    rf_num_trees = 5L,
    project_root = ".",
    n_cores = 1L
) {
  if (!requireNamespace("batchtools", quietly = TRUE)) {
    stop("Package 'batchtools' is required for the socket smoke test.", call. = FALSE)
  }
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for the socket smoke test.", call. = FALSE)
  }

  td <- tempfile(pattern = "signal_sanity_bt_smoke_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)

  bootstrap <- signal_sanity_write_worker_bootstrap(
    destination = file.path(td, "worker_bootstrap.R"),
    project_root = project_root
  )

  reg <- batchtools::makeRegistry(
    file.dir = file.path(td, "registry"),
    packages = "ranger",
    source = bootstrap,
    seed = 1L,
    make.default = FALSE
  )

  batchtools::batchMap(
    fun = signal_sanity_bt_worker,
    task_row = list(as.data.frame(task_row, stringsAsFactors = FALSE)),
    more.args = list(
      settings = settings,
      design = design,
      rf_num_trees = as.integer(rf_num_trees),
      analysis_n_test = as.integer(analysis_n_test)
    ),
    reg = reg
  )

  reg <- set_batchtools_cluster_function(
    reg = reg,
    n_cores = as.integer(n_cores),
    cluster_function = "socket"
  )

  batchtools::submitJobs(reg = reg)
  batchtools::waitForJobs(reg = reg)
  try(batchtools::syncRegistry(reg = reg), silent = TRUE)

  err <- signal_sanity_job_ids(batchtools::findErrors(reg = reg))
  not_done <- signal_sanity_job_ids(batchtools::findNotDone(reg = reg))
  if (length(err) > 0L || length(not_done) > 0L) {
    stop(
      "Batchtools/socket smoke test failed before the full run. errors=",
      length(err), ", not_done=", length(not_done),
      ". Inspect the smoke-test registry messages, the 9,600-job run was NOT started.",
      call. = FALSE
    )
  }

  ans <- batchtools::reduceResultsList(reg = reg, missing.val = NULL)
  ans <- signal_sanity_bind_results(ans)
  if (nrow(ans) != 1L || !isTRUE(ans$ok[[1L]])) {
    msg <- if (nrow(ans) == 1L && "error_message" %in% names(ans)) ans$error_message[[1L]] else "unknown"
    stop("Batchtools/socket smoke test returned an invalid scientific result: ", msg, call. = FALSE)
  }

  return(invisible(ans))
}

# Run or resume the full-grid signal check with 50 repetitions and RF with 200 trees.
# Inputs: pilot size, analysis size, worker/seed options and output/registry paths.
# Returns the registry, tasks, results and descriptive summaries invisibly.
# save_outputs writes tables, a summary object and session information.
# Failure rows are retained before stopping. This is diagnostic, not a hypothesis test.
run_signal_sanity_full_grid <- function(
    R = 50L,
    n_test = 5000L,
    rf_num_trees = 200L,
    n_cores = 8L,
    seed_base = 1300000000L,
    output_dir = "outputs/signal_sanity",
    registry_dir = file.path(output_dir, "registry_R50_corrected"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = TRUE,
    wait = TRUE,
    save_outputs = TRUE
) {

  signal_sanity_assert_core()
  if (!requireNamespace("batchtools", quietly = TRUE)) stop("Package 'batchtools' is required.", call. = FALSE)
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' is required.", call. = FALSE)

  R <- as.integer(R)
  n_test <- as.integer(n_test)
  rf_num_trees <- as.integer(rf_num_trees)
  n_cores <- as.integer(n_cores)
  seed_base <- as.integer(seed_base)
  if (R != 50L) warning("This reviewed diagnostic was planned for R=50, requested R=", R, ".")
  if (rf_num_trees != 200L) stop("Full diagnostic must use rf_num_trees = 200 to match the precision pilot.", call. = FALSE)
  if (n_cores < 1L) stop("n_cores must be >= 1.", call. = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  design <- signal_sanity_prepare_design(n_test = n_test)
  settings <- signal_sanity_make_settings(design, expected_n_settings = 192L)
  tasks <- signal_sanity_make_tasks(settings, R = R, seed_base = seed_base)

  signature_source_files <- c("R/pilot_core.R", "R/precision_pilot.R", "R/signal_sanity.R")
  missing_source <- signature_source_files[!file.exists(signature_source_files)]
  if (length(missing_source) > 0L) stop("Missing source file(s): ", paste(missing_source, collapse = ", "), call. = FALSE)

  project_root_abs <- normalizePath(".", winslash = "/", mustWork = TRUE)
  worker_bootstrap <- signal_sanity_write_worker_bootstrap(
    destination = file.path(output_dir, "signal_sanity_worker_bootstrap.R"),
    project_root = project_root_abs
  )

  signature <- signal_sanity_signature(
    R = R, n_test = n_test, rf_num_trees = rf_num_trees, n_cores = n_cores,
    seed_base = seed_base, settings = settings, source_files = signature_source_files
  )

  registry_exists <- dir.exists(registry_dir) && file.exists(file.path(registry_dir, "registry.rds"))
  signature_file <- file.path(registry_dir, "signal_sanity_signature.rds")

  if (registry_exists && overwrite_registry) {
    unlink(registry_dir, recursive = TRUE, force = TRUE)
    registry_exists <- FALSE
  }

  if (registry_exists) {
    if (!resume) stop("Registry already exists and resume=FALSE: ", registry_dir, call. = FALSE)
    if (!file.exists(signature_file)) stop("Existing registry has no signature file, refusing unsafe resume.", call. = FALSE)
    old_signature <- readRDS(signature_file)
    if (!identical(old_signature, signature)) {
      stop("Existing registry signature does not match this run. Use a new registry path, do not overwrite completed scientific results.", call. = FALSE)
    }
    reg <- batchtools::loadRegistry(file.dir = registry_dir, writeable = TRUE, make.default = FALSE)
    try(batchtools::syncRegistry(reg = reg), silent = TRUE)
  } else {
    dir.create(dirname(registry_dir), recursive = TRUE, showWarnings = FALSE)
    reg <- batchtools::makeRegistry(
      file.dir = registry_dir,
      packages = "ranger",
      source = worker_bootstrap,
      seed = 1L,
      make.default = FALSE
    )
    task_list <- split(tasks, seq_len(nrow(tasks)))
    batchtools::batchMap(
      fun = signal_sanity_bt_worker,
      task_row = task_list,
      more.args = list(settings = settings, design = design, rf_num_trees = rf_num_trees, analysis_n_test = n_test),
      reg = reg
    )
    dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(signature, signature_file, version = 2)
    saveRDS(settings, file.path(registry_dir, "signal_sanity_settings.rds"), version = 2)
    saveRDS(tasks, file.path(registry_dir, "signal_sanity_tasks.rds"), version = 2)
  }

  reg <- set_batchtools_cluster_function(
    reg = reg,
    n_cores = n_cores,
    cluster_function = "socket"
  )

  on_system <- signal_sanity_job_ids(batchtools::findOnSystem(reg = reg))
  if (length(on_system) > 0L) {
    stop(
      "Registry still reports ", length(on_system), " job(s) as 'on system'. Do not start a second copy. ",
      "After a confirmed server restart these stale jobs can be reset safely before resuming.",
      call. = FALSE
    )
  }

  err_ids <- signal_sanity_job_ids(batchtools::findErrors(reg = reg))
  if (length(err_ids) > 0L && isTRUE(resubmit_errors)) {
    batchtools::resetJobs(ids = err_ids, reg = reg)
  }

  todo <- signal_sanity_job_ids(batchtools::findNotDone(reg = reg))
  if (length(todo) > 0L) {
    message("Submitting ", length(todo), " unfinished signal-sanity jobs out of ", nrow(tasks), ".")
    batchtools::submitJobs(ids = todo, reg = reg)
    if (isTRUE(wait)) batchtools::waitForJobs(ids = todo, reg = reg)
  }

  if (!isTRUE(wait)) {
    return(invisible(list(registry = reg, settings = settings, tasks = tasks, signature = signature)))
  }

  try(batchtools::syncRegistry(reg = reg), silent = TRUE)

  err_after <- signal_sanity_job_ids(batchtools::findErrors(reg = reg))
  if (length(err_after) > 0L && isTRUE(resubmit_errors)) {
    message("Retrying ", length(err_after), " batchtools error job(s) once.")
    batchtools::resetJobs(ids = err_after, reg = reg)
    batchtools::submitJobs(ids = err_after, reg = reg)
    batchtools::waitForJobs(ids = err_after, reg = reg)
    try(batchtools::syncRegistry(reg = reg), silent = TRUE)
  }

  remaining <- signal_sanity_job_ids(batchtools::findNotDone(reg = reg))
  if (length(remaining) > 0L) {
    stop(length(remaining), " jobs are still not done after wait/retry. Inspect the registry before resuming.", call. = FALSE)
  }

  result_list <- batchtools::reduceResultsList(reg = reg, missing.val = NULL)
  results <- signal_sanity_bind_results(result_list)

  failure_records <- results[is.na(results$ok) | !results$ok, , drop = FALSE]
  if (nrow(failure_records) > 0L && isTRUE(save_outputs)) {
    utils::write.csv(
      failure_records, file.path(output_dir, "signal_sanity_failures.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      results, file.path(output_dir, "signal_sanity_task_results_partial.csv"),
      row.names = FALSE
    )
  }
  if (nrow(results) != nrow(tasks)) {
    stop("Collected ", nrow(results), " results but expected ", nrow(tasks), ".", call. = FALSE)
  }
  results <- results[order(results$task_id), , drop = FALSE]
  expected_ids <- seq_len(nrow(tasks))
  if (!identical(as.integer(results$task_id), as.integer(expected_ids))) {
    stop("Task-result IDs are incomplete or out of alignment.", call. = FALSE)
  }

  failed <- results[is.na(results$ok) | !results$ok, , drop = FALSE]
  successful <- results[!is.na(results$ok) & results$ok, , drop = FALSE]

  if (nrow(successful) > 0L) {
    setting_summary <- signal_sanity_setting_summary(successful)
    by_snr <- signal_sanity_by_snr(setting_summary)
    overview <- signal_sanity_overview(results, setting_summary, expected_jobs = nrow(tasks))
    low_gain <- setting_summary[order(setting_summary$mean_rf_gain_vs_intercept), , drop = FALSE]
  } else {
    setting_summary <- data.frame()
    by_snr <- data.frame()
    overview <- data.frame()
    low_gain <- data.frame()
  }

  if (isTRUE(save_outputs)) {
    utils::write.csv(results, file.path(output_dir, "signal_sanity_task_results.csv"), row.names = FALSE)
    utils::write.csv(settings, file.path(output_dir, "signal_sanity_design_snapshot.csv"), row.names = FALSE)
    utils::write.csv(setting_summary, file.path(output_dir, "signal_sanity_setting_summary.csv"), row.names = FALSE)
    utils::write.csv(by_snr, file.path(output_dir, "signal_sanity_by_snr.csv"), row.names = FALSE)
    utils::write.csv(overview, file.path(output_dir, "signal_sanity_overview.csv"), row.names = FALSE)
    utils::write.csv(utils::head(low_gain, 30L), file.path(output_dir, "signal_sanity_lowest_rf_gain_settings.csv"), row.names = FALSE)
    utils::write.csv(failed, file.path(output_dir, "signal_sanity_failures.csv"), row.names = FALSE)
    saveRDS(
      list(
        version = signal_sanity_version,
        signature = signature,
        settings = settings,
        tasks = tasks,
        results = results,
        setting_summary = setting_summary,
        by_snr = by_snr,
        overview = overview,
        registry_status = batchtools::getStatus(reg = reg)
      ),
      file.path(output_dir, "signal_sanity_summary.rds"),
      version = 2
    )
    capture.output(sessionInfo(), file = file.path(output_dir, "signal_sanity_sessionInfo.txt"))
  }

  if (nrow(failed) > 0L) {
    stop(
      "The batchtools run finished, but ", nrow(failed), " task(s) returned ok=FALSE. ",
      "Outputs were saved, inspect signal_sanity_failures.csv before interpreting results.",
      call. = FALSE
    )
  }

  return(invisible(list(
    registry = reg,
    settings = settings,
    tasks = tasks,
    results = results,
    setting_summary = setting_summary,
    by_snr = by_snr,
    overview = overview,
    signature = signature
  )))
}
