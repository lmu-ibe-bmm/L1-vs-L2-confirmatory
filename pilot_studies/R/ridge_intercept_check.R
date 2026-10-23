# Ridge versus the training-mean baseline in the 12 n = 1000, p = 5000, rho = 0 settings.

ridge_intercept_check_version <- "1.0.0"

ridge_intercept_required_scientific_functions <- c(
  "make_design", "make_all_settings", "make_beta", "simulate_data",
  "make_cv_foldid", "fit_glmnet"
)

# Check that the shared DGM and fitting functions are loaded.
# Inputs: none. Returns invisible TRUE, or stops if a dependency is missing.
ridge_intercept_assert_scientific_core <- function() {
  missing <- ridge_intercept_required_scientific_functions[
    !vapply(
      ridge_intercept_required_scientific_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
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

# Check that the socket-worker setup is available in the calling session.
# Inputs: none. Returns invisible TRUE or stops before creating jobs.
ridge_intercept_assert_master_infrastructure <- function() {
  if (!exists("set_batchtools_cluster_function", mode = "function", inherits = TRUE)) {
    stop(
      "set_batchtools_cluster_function() is missing. ",
      "Source R/precision_pilot.R in the master R session.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}

# Validate one task field before converting its type.
# Inputs: x and its field name. Returns the single non-missing value.
ridge_intercept_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x)) {
    stop("Expected one non-missing value for '", name, "'.", call. = FALSE)
  }
  return(x[[1L]])
}

# Use the shared Gaussian DGM and glmnet settings for this diagnostic.
# Inputs: none. Returns the checked design with positive noise variance and valid folds.
ridge_intercept_prepare_design <- function() {
  ridge_intercept_assert_scientific_core()

  design <- make_design()

  if (!is.list(design)) {
    stop("make_design() did not return a list-like design object.", call. = FALSE)
  }

  if (is.null(design$sigma2) || length(design$sigma2) != 1L ||
    !is.finite(design$sigma2) || design$sigma2 <= 0) {
    stop("design$sigma2 must exist and be one positive finite value.", call. = FALSE)
  }

  if (is.null(design$glmnet_nfolds) || length(design$glmnet_nfolds) != 1L ||
    is.na(design$glmnet_nfolds) || as.integer(design$glmnet_nfolds) < 3L) {
    stop("design$glmnet_nfolds must exist and be an integer >= 3.", call. = FALSE)
  }

  return(design)
}

# Build and validate the full design before selecting the diagnostic subset.
# Inputs: design and expected_n_settings. Returns settings ordered by their original IDs.
ridge_intercept_full_settings <- function(design, expected_n_settings = 192L) {
  settings <- as.data.frame(make_all_settings(design), stringsAsFactors = FALSE)

  if (nrow(settings) != as.integer(expected_n_settings)) {
    stop(
      "Expected ", expected_n_settings,
      " settings in the current Registered-Report grid, but make_all_settings() returned ",
      nrow(settings),
      ". Stopping rather than silently testing a different design.",
      call. = FALSE
    )
  }

  required <- c("n", "p_over_n", "p", "rho", "a", "k", "snr", "b")
  missing <- setdiff(required, names(settings))
  if (length(missing) > 0L) {
    stop(
      "Settings are missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"setting_id" %in% names(settings)) {
    settings$setting_id <- seq_len(nrow(settings))
  }
  settings$setting_id <- as.integer(settings$setting_id)

  if (anyNA(settings$setting_id) || anyDuplicated(settings$setting_id)) {
    stop("setting_id must be complete and unique.", call. = FALSE)
  }

  if (!identical(sort(settings$setting_id), seq_len(nrow(settings)))) {
    stop("Expected setting_id to be exactly 1,...,192.", call. = FALSE)
  }

  expected_p <- as.integer(round(settings$n * settings$p_over_n))
  if (any(expected_p != as.integer(settings$p))) {
    bad <- which(expected_p != as.integer(settings$p))[1L]
    stop(
      "Inconsistent n, p/n and p in setting_id ",
      settings$setting_id[bad],
      ".",
      call. = FALSE
    )
  }

  return(settings[order(settings$setting_id), , drop = FALSE])
}

# Select IDs 145 to 156: n = 1000, p = 5000 and rho = 0.
# Inputs: full_settings. Returns all 12 a, k and SNR combinations in that subset.
ridge_intercept_select_settings <- function(full_settings) {

  keep <- as.integer(full_settings$p) == 5000L &
    abs(as.numeric(full_settings$rho)) < 1e-12

  selected <- full_settings[keep, , drop = FALSE]
  selected <- selected[order(selected$setting_id), , drop = FALSE]

  if (nrow(selected) != 12L) {
    stop(
      "Expected exactly 12 p=5000, rho=0 settings, but found ",
      nrow(selected),
      ". The current study grid differs from the reviewed design.",
      call. = FALSE
    )
  }

  if (!identical(sort(unique(as.integer(selected$n))), 1000L)) {
    stop("Selected settings are expected to have n=1000.", call. = FALSE)
  }

  if (!isTRUE(all(abs(as.numeric(selected$p_over_n) - 5) < 1e-12))) {
    stop("Selected settings are expected to have p/n=5.", call. = FALSE)
  }

  if (!identical(sort(unique(as.numeric(selected$a))), c(0.25, 0.75))) {
    stop("Expected a in {0.25, 0.75} for selected settings.", call. = FALSE)
  }

  if (!identical(sort(unique(as.integer(selected$k))), c(1L, 5L, 10L))) {
    stop("Expected k in {1, 5, 10} for selected settings.", call. = FALSE)
  }

  if (!identical(sort(unique(as.numeric(selected$snr))), c(0.5, 2))) {
    stop("Expected SNR in {0.5, 2} for selected settings.", call. = FALSE)
  }

  expected_combo <- expand.grid(
    a = c(0.25, 0.75),
    k = c(1L, 5L, 10L),
    snr = c(0.5, 2),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  observed_key <- paste(selected$a, selected$k, selected$snr, sep = "|")
  expected_key <- paste(expected_combo$a, expected_combo$k, expected_combo$snr, sep = "|")

  if (!setequal(observed_key, expected_key)) {
    stop(
      "The selected settings do not contain the expected full a x k x SNR cross.",
      call. = FALSE
    )
  }

  if (!identical(as.integer(selected$setting_id), 145:156)) {
    stop(
      "The correct design values are present, but setting IDs are not 145:156. ",
      "Stopping to prevent accidental comparison of differently ordered grids.",
      call. = FALSE
    )
  }

  return(selected)
}

# Create one reproducible task for each selected setting and repetition.
# Inputs: selected_settings, R and seed_base. Returns task IDs and corrected component-seed bases.
ridge_intercept_make_tasks <- function(
    selected_settings,
    R = 50L,
    seed_base = 1700000000L
) {

  if (length(R) != 1L || !is.finite(R) || R != floor(R) ||
    R < 1 || R > 9999) {
    stop("R must be one integer from 1 to 9999.", call. = FALSE)
  }
  if (nrow(selected_settings) < 1L || anyDuplicated(selected_settings$setting_id)) {
    stop("Settings must have distinct setting IDs.", call. = FALSE)
  }
  make_task_seeds(selected_settings$setting_id, rep(R, nrow(selected_settings)), seed_base)
  R <- as.integer(R)

  grid <- expand.grid(
    setting_id = as.integer(selected_settings$setting_id),
    run_id = seq_len(R),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid <- grid[order(grid$setting_id, grid$run_id), , drop = FALSE]
  rownames(grid) <- NULL
  grid$task_id <- seq_len(nrow(grid))

  grid$seed <- make_task_seeds(grid$setting_id, grid$run_id, seed_base)

  component_seeds <- unlist(

    lapply(1:4, function(off) return(as.double(grid$seed) + off)),
    use.names = FALSE
  )

  if (anyDuplicated(component_seeds)) {
    stop("Internal error: component-seed collision detected.", call. = FALSE)
  }

  return(grid[, c("task_id", "setting_id", "run_id", "seed")])
}

# Check the population signal variance under block equicorrelation.
# Inputs: beta, rho and block size b. Returns beta transpose Sigma beta.
ridge_intercept_signal_variance <- function(beta, rho, b) {
  beta <- as.numeric(beta)
  rho <- as.numeric(rho)
  b <- as.integer(b)

  if (length(beta) < 1L || any(!is.finite(beta))) {
    stop("beta must be finite and non-empty.", call. = FALSE)
  }

  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
    stop("rho must be in [0,1).", call. = FALSE)
  }

  if (length(b) != 1L || is.na(b) || b < 1L || length(beta) %% b != 0L) {
    stop("b must be a positive divisor of length(beta).", call. = FALSE)
  }

  block_id <- rep(seq_len(length(beta) / b), each = b)
  block_sums <- vapply(split(beta, block_id), sum, numeric(1))

  return((1 - rho) * sum(beta^2) + rho * sum(block_sums^2))
}

# Calculate test MSE from aligned outcomes y and predictions pred.
# Returns one mean squared prediction error.
ridge_intercept_mse <- function(y, pred) {
  return(mean((as.numeric(y) - as.numeric(pred))^2))
}

# Read a finite scalar from the shared glmnet fit result.
# Inputs: fit, candidate field names and label. Returns the value or stops.
ridge_intercept_extract_fit_value <- function(fit, candidates, label) {
  for (nm in candidates) {
    if (!is.null(fit[[nm]])) {
      x <- fit[[nm]]
      if (length(x) == 1L && !is.na(x) && is.finite(as.numeric(x))) {
        return(as.numeric(x))
      }
    }
  }

  stop(
    "fit_glmnet() result did not contain a finite scalar for ",
    label,
    ". Tried: ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

# Check Ridge against the training-mean baseline on the same test sample.
# Inputs: task_row, full/selected settings, design and analysis_n_test.
# Returns one row with MSE, relative gains, SNR and task identifiers.
# Gain = 1 - model MSE / baseline MSE. Positive values indicate improvement.
# The known-coefficient prediction is a signal reference. Failures return ok = FALSE.
ridge_intercept_one_task <- function(
    task_row,
    full_settings,
    selected_settings,
    design,
    analysis_n_test = 5000L
) {
  started <- proc.time()[[3L]]

  task_id <- setting_id <- run_id <- seed <- NA_integer_
  analysis_n_test <- as.integer(analysis_n_test)

  if (length(analysis_n_test) != 1L || is.na(analysis_n_test) ||
    analysis_n_test < 100L) {
    stop("analysis_n_test must be one integer >= 100.", call. = FALSE)
  }

  return(tryCatch({

    task_row <- as.data.frame(task_row, stringsAsFactors = FALSE)

    if (nrow(task_row) != 1L) {
      stop("task_row must contain exactly one row.")
    }

    task_id <- as.integer(
      ridge_intercept_scalar(task_row$task_id, "task_id")
    )
    setting_id <- as.integer(
      ridge_intercept_scalar(task_row$setting_id, "setting_id")
    )
    run_id <- as.integer(
      ridge_intercept_scalar(task_row$run_id, "run_id")
    )
    seed <- as.integer(
      ridge_intercept_scalar(task_row$seed, "seed")
    )

    idx_full <- match(setting_id, full_settings$setting_id)
    idx_sel <- match(setting_id, selected_settings$setting_id)

    if (is.na(idx_full) || is.na(idx_sel)) {
      stop("Unknown or non-selected setting_id: ", setting_id)
    }

    setting <- full_settings[idx_full, , drop = FALSE]

    worker_row <- setting
    worker_row$task_id <- task_id
    worker_row$run_id <- run_id
    worker_row$seed <- seed

    beta_obj <- make_beta(
      setting_row = worker_row,
      seed = seed + 1L
    )

    if (is.null(beta_obj$beta)) {
      stop("make_beta() result has no $beta component.")
    }

    beta <- as.numeric(beta_obj$beta)

    dat_native <- simulate_data(
      setting_row = worker_row,
      beta_object = beta_obj,
      seed = seed + 2L,
      design = design
    )

    required_dat <- c("X_train", "y_train", "X_test", "y_test")
    missing_dat <- setdiff(required_dat, names(dat_native))

    if (length(missing_dat) > 0L) {
      stop(
        "simulate_data() result is missing: ",
        paste(missing_dat, collapse = ", ")
      )
    }

    X_train <- dat_native$X_train
    y_train <- as.numeric(dat_native$y_train)
    X_test_native <- dat_native$X_test
    y_test_native <- as.numeric(dat_native$y_test)

    p <- ncol(X_train)

    if (is.null(p) || p != length(beta) || ncol(X_test_native) != p) {
      stop("X/beta dimension mismatch.")
    }

    if (length(y_train) != nrow(X_train) ||
      length(y_test_native) != nrow(X_test_native)) {
      stop("X/y dimension mismatch.")
    }

    native_n_test <- nrow(X_test_native)

    if (native_n_test < analysis_n_test) {
      stop(
        "Project-native test sample has only ",
        native_n_test,
        " rows, but this diagnostic requests ",
        analysis_n_test,
        "."
      )
    }

    analysis_idx <- seq_len(analysis_n_test)
    X_test <- X_test_native[analysis_idx, , drop = FALSE]
    y_test <- y_test_native[analysis_idx]

    dat <- dat_native
    dat$X_test <- X_test
    dat$y_test <- y_test

    nfolds <- as.integer(design$glmnet_nfolds)
    fold_seed <- seed + 3L

    foldid <- make_cv_foldid(
      n_train = nrow(X_train),
      nfolds = nfolds,
      seed = fold_seed
    )

    if (length(foldid) != nrow(X_train) || anyNA(foldid)) {
      stop("make_cv_foldid() returned invalid fold IDs.")
    }

    fold_levels <- sort(unique(as.integer(foldid)))

    if (!identical(fold_levels, seq_len(nfolds))) {
      stop(
        "CV fold IDs are not exactly 1,...,glmnet_nfolds. Observed: ",
        paste(fold_levels, collapse = ",")
      )
    }

    pred_intercept <- rep(mean(y_train), length(y_test))
    pred_oracle <- drop(X_test %*% beta)

    mse_intercept <- ridge_intercept_mse(
      y_test,
      pred_intercept
    )
    mse_oracle <- ridge_intercept_mse(
      y_test,
      pred_oracle
    )

    if (!is.finite(mse_intercept) || mse_intercept <= 0 ||
      !is.finite(mse_oracle) || mse_oracle <= 0) {
      stop("Non-finite or non-positive baseline MSE encountered.")
    }

    ridge_fit <- fit_glmnet(
      dat = dat,
      alpha = 0,
      foldid = foldid,
      design = design
    )

    mse_ridge <- ridge_intercept_extract_fit_value(
      ridge_fit,
      "mse",
      "Ridge test MSE"
    )
    lambda_min <- ridge_intercept_extract_fit_value(
      ridge_fit,
      c("lambda_min", "lambda.min"),
      "Ridge lambda.min"
    )
    lambda_1se <- ridge_intercept_extract_fit_value(
      ridge_fit,
      c("lambda_1se", "lambda.1se"),
      "Ridge lambda.1se"
    )

    if (!is.finite(mse_ridge) || mse_ridge <= 0) {
      stop("Ridge returned non-finite or non-positive test MSE.")
    }

    target_snr <- as.numeric(setting$snr[[1L]])
    rho <- as.numeric(setting$rho[[1L]])
    b <- as.integer(setting$b[[1L]])
    sigma2 <- as.numeric(design$sigma2)

    signal_var_analytic <- ridge_intercept_signal_variance(
      beta,
      rho = rho,
      b = b
    )
    analytic_snr <- signal_var_analytic / sigma2
    analytic_snr_rel_error <- analytic_snr / target_snr - 1

    eta_test <- pred_oracle
    eps_test <- y_test - eta_test
    empirical_snr_test <- stats::var(eta_test) / stats::var(eps_test)

    ridge_gain <- 1 - mse_ridge / mse_intercept
    ridge_difference <- mse_intercept - mse_ridge
    oracle_gain <- 1 - mse_oracle / mse_intercept
    oracle_denom <- mse_intercept - mse_oracle

    ridge_fraction_oracle_gain <-
      if (is.finite(oracle_denom) && oracle_denom > 0) {
        (mse_intercept - mse_ridge) / oracle_denom
      } else {
        NA_real_
      }

    data.frame(
      ok = TRUE,
      error_message = NA_character_,
      task_id = task_id,
      setting_id = setting_id,
      run_id = run_id,
      seed = seed,
      beta_seed = seed + 1L,
      data_seed = seed + 2L,
      cv_fold_seed = fold_seed,
      model = "Ridge_alpha_0",
      alpha = 0,
      n = as.integer(setting$n[[1L]]),
      p_over_n = as.numeric(setting$p_over_n[[1L]]),
      p = as.integer(setting$p[[1L]]),
      rho = rho,
      a = as.numeric(setting$a[[1L]]),
      k = as.integer(setting$k[[1L]]),
      snr = target_snr,
      b = b,
      n_test = as.integer(length(y_test)),
      n_test_generated = as.integer(native_n_test),
      glmnet_nfolds = nfolds,
      cv_fold_min_size = min(
        tabulate(as.integer(foldid), nbins = nfolds)
      ),
      cv_fold_max_size = max(
        tabulate(as.integer(foldid), nbins = nfolds)
      ),
      analytic_snr = analytic_snr,
      analytic_snr_rel_error = analytic_snr_rel_error,
      empirical_snr_test = empirical_snr_test,
      population_oracle_gain_reference =
        target_snr / (1 + target_snr),
      mse_intercept = mse_intercept,
      mse_oracle = mse_oracle,
      mse_ridge = mse_ridge,
      ridge_mse_difference_vs_intercept = ridge_difference,
      ridge_gain_vs_intercept = ridge_gain,
      oracle_gain_vs_intercept = oracle_gain,
      ridge_fraction_oracle_gain = ridge_fraction_oracle_gain,
      ridge_beats_intercept = isTRUE(mse_ridge < mse_intercept),
      oracle_beats_intercept = isTRUE(mse_oracle < mse_intercept),
      oracle_mse_over_sigma2 = mse_oracle / sigma2,
      intercept_mse_over_population_reference =
        mse_intercept / (sigma2 * (1 + target_snr)),
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      elapsed_sec = proc.time()[[3L]] - started,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    return(data.frame(
      ok = FALSE,
      error_message = conditionMessage(e),
      task_id =
        task_id,
      setting_id =
        setting_id,
      run_id =
        run_id,
      seed =
        seed,
      model = "Ridge_alpha_0",
      alpha = 0,
      elapsed_sec = proc.time()[[3L]] - started,
      stringsAsFactors = FALSE
    ))
  }))
}

# Batchtools entry point with the same inputs and result as the single-task function.
# Delegates the diagnostic calculation without changing the DGM or fit options.
ridge_intercept_bt_worker <- function(
    task_row,
    full_settings,
    selected_settings,
    design,
    analysis_n_test = 5000L
) {
  return(ridge_intercept_one_task(
    task_row = task_row,
    full_settings = full_settings,
    selected_settings = selected_settings,
    design = design,
    analysis_n_test = analysis_n_test
  ))
}

# Summarize finite repetition-level values with a descriptive t interval.
# Inputs: x and confidence level. Returns mean, Monte Carlo SE and bounds.
# The SE is sd(x) / sqrt(n). Fewer than two values give missing bounds.
ridge_intercept_mean_ci <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(c(
      mean = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  }

  m <- mean(x)

  if (n == 1L) {
    return(c(
      mean = m,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  }

  se <- stats::sd(x) / sqrt(n)
  crit <- stats::qt((1 + level) / 2, df = n - 1L)

  return(c(
    mean = m,
    se = se,
    lower = m - crit * se,
    upper = m + crit * se
  ))
}

# Aggregate successful repetitions within each setting.
# Inputs: results. Returns MSE/gain summaries, SNR checks and descriptive t intervals.
# These are summaries of repetition-level gains, not the confirmatory ratio of means.
ridge_intercept_setting_summary <- function(results) {
  ok <- results[!is.na(results$ok) & results$ok, , drop = FALSE]

  if (nrow(ok) == 0L) {
    stop("No successful task results to summarise.", call. = FALSE)
  }

  groups <- split(ok, ok$setting_id)

  rows <- lapply(groups, function(d) {
    gain_ci <- ridge_intercept_mean_ci(
      d$ridge_gain_vs_intercept
    )
    diff_ci <- ridge_intercept_mean_ci(
      d$ridge_mse_difference_vs_intercept
    )
    fraction_ci <- ridge_intercept_mean_ci(
      d$ridge_fraction_oracle_gain
    )
    snr_ci <- ridge_intercept_mean_ci(
      d$empirical_snr_test
    )

    return(data.frame(
      setting_id = as.integer(d$setting_id[1L]),
      n_success = nrow(d),
      model = as.character(d$model[1L]),
      alpha = as.numeric(d$alpha[1L]),
      n = as.integer(d$n[1L]),
      p_over_n = as.numeric(d$p_over_n[1L]),
      p = as.integer(d$p[1L]),
      rho = as.numeric(d$rho[1L]),
      a = as.numeric(d$a[1L]),
      k = as.integer(d$k[1L]),
      snr = as.numeric(d$snr[1L]),
      n_test = as.integer(d$n_test[1L]),
      min_n_test_generated =
        min(as.integer(d$n_test_generated)),
      max_n_test_generated =
        max(as.integer(d$n_test_generated)),
      glmnet_nfolds = as.integer(d$glmnet_nfolds[1L]),
      min_cv_fold_size =
        min(as.integer(d$cv_fold_min_size)),
      max_cv_fold_size =
        max(as.integer(d$cv_fold_max_size)),
      mean_analytic_snr =
        mean(d$analytic_snr),
      max_abs_analytic_snr_rel_error =
        max(abs(d$analytic_snr_rel_error)),
      mean_empirical_snr_test =
        snr_ci[["mean"]],
      empirical_snr_test_ci_lower =
        snr_ci[["lower"]],
      empirical_snr_test_ci_upper =
        snr_ci[["upper"]],
      mean_mse_intercept =
        mean(d$mse_intercept),
      sd_mse_intercept =
        stats::sd(d$mse_intercept),
      mean_mse_oracle =
        mean(d$mse_oracle),
      sd_mse_oracle =
        stats::sd(d$mse_oracle),
      mean_mse_ridge =
        mean(d$mse_ridge),
      sd_mse_ridge =
        stats::sd(d$mse_ridge),
      mean_ridge_mse_difference_vs_intercept =
        diff_ci[["mean"]],
      ridge_mse_difference_ci_lower =
        diff_ci[["lower"]],
      ridge_mse_difference_ci_upper =
        diff_ci[["upper"]],
      mean_ridge_gain_vs_intercept =
        gain_ci[["mean"]],
      median_ridge_gain_vs_intercept =
        stats::median(d$ridge_gain_vs_intercept),
      sd_ridge_gain_vs_intercept =
        stats::sd(d$ridge_gain_vs_intercept),
      ridge_gain_ci_lower =
        gain_ci[["lower"]],
      ridge_gain_ci_upper =
        gain_ci[["upper"]],
      q10_ridge_gain_vs_intercept =
        as.numeric(
          stats::quantile(
            d$ridge_gain_vs_intercept,
            0.10,
            names = FALSE
          )
        ),
      q90_ridge_gain_vs_intercept =
        as.numeric(
          stats::quantile(
            d$ridge_gain_vs_intercept,
            0.90,
            names = FALSE
          )
        ),
      prop_ridge_beats_intercept =
        mean(d$ridge_beats_intercept),
      mean_oracle_gain_vs_intercept =
        mean(d$oracle_gain_vs_intercept),
      mean_population_oracle_gain_reference =
        mean(d$population_oracle_gain_reference),
      mean_ridge_fraction_oracle_gain =
        fraction_ci[["mean"]],
      median_ridge_fraction_oracle_gain =
        stats::median(
          d$ridge_fraction_oracle_gain,
          na.rm = TRUE
        ),
      ridge_fraction_oracle_gain_ci_lower =
        fraction_ci[["lower"]],
      ridge_fraction_oracle_gain_ci_upper =
        fraction_ci[["upper"]],
      median_lambda_min =
        stats::median(d$lambda_min),
      median_lambda_1se =
        stats::median(d$lambda_1se),
      mean_oracle_mse_over_sigma2 =
        mean(d$oracle_mse_over_sigma2),
      mean_intercept_mse_over_population_reference =
        mean(d$intercept_mse_over_population_reference),
      mean_elapsed_sec =
        mean(d$elapsed_sec),
      stringsAsFactors = FALSE
    ))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  return(out[order(out$setting_id), , drop = FALSE])
}

# Describe heterogeneity across settings within each SNR level.
# Inputs: setting_summary. Returns setting counts and the distribution of estimated gains.
ridge_intercept_by_snr <- function(setting_summary) {
  groups <- split(
    setting_summary,
    setting_summary$snr
  )

  rows <- lapply(groups, function(d) {
    return(data.frame(
      snr = as.numeric(d$snr[1L]),
      n_settings = nrow(d),
      min_mean_ridge_gain_vs_intercept =
        min(d$mean_ridge_gain_vs_intercept),
      q25_mean_ridge_gain_vs_intercept =
        as.numeric(
          stats::quantile(
            d$mean_ridge_gain_vs_intercept,
            0.25,
            names = FALSE
          )
        ),
      median_mean_ridge_gain_vs_intercept =
        stats::median(d$mean_ridge_gain_vs_intercept),
      q75_mean_ridge_gain_vs_intercept =
        as.numeric(
          stats::quantile(
            d$mean_ridge_gain_vs_intercept,
            0.75,
            names = FALSE
          )
        ),
      max_mean_ridge_gain_vs_intercept =
        max(d$mean_ridge_gain_vs_intercept),
      min_prop_ridge_beats_intercept =
        min(d$prop_ridge_beats_intercept),
      n_settings_ridge_gain_ci_lower_gt_0 =
        sum(d$ridge_gain_ci_lower > 0, na.rm = TRUE),
      n_settings_mean_ridge_gain_le_0 =
        sum(d$mean_ridge_gain_vs_intercept <= 0, na.rm = TRUE),
      min_mean_ridge_fraction_oracle_gain =
        min(d$mean_ridge_fraction_oracle_gain, na.rm = TRUE),
      median_mean_ridge_fraction_oracle_gain =
        stats::median(
          d$mean_ridge_fraction_oracle_gain,
          na.rm = TRUE
        ),
      max_abs_analytic_snr_rel_error =
        max(d$max_abs_analytic_snr_rel_error),
      stringsAsFactors = FALSE
    ))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  return(out[order(out$snr), , drop = FALSE])
}

# Describe run completeness, gain ranges and the least improved setting.
# Inputs: results, setting_summary and expected_jobs. Returns one diagnostic overview.
ridge_intercept_overview <- function(
    results,
    setting_summary,
    expected_jobs
) {
  n_success <- sum(!is.na(results$ok) & results$ok)
  n_failed <- nrow(results) - n_success

  worst <- setting_summary[
    which.min(setting_summary$mean_ridge_gain_vs_intercept),
    ,
    drop = FALSE
  ]

  return(data.frame(
    version = ridge_intercept_check_version,
    expected_jobs = as.integer(expected_jobs),
    n_results = nrow(results),
    n_success = n_success,
    n_failed = n_failed,
    n_settings = nrow(setting_summary),
    repetitions_per_setting =
      if (nrow(setting_summary) > 0L) {
        min(setting_summary$n_success)
      } else {
        NA_integer_
      },
    n_settings_mean_ridge_gain_gt_0 =
      sum(setting_summary$mean_ridge_gain_vs_intercept > 0, na.rm = TRUE),
    n_settings_ridge_gain_ci_lower_gt_0 =
      sum(setting_summary$ridge_gain_ci_lower > 0, na.rm = TRUE),
    n_settings_mean_ridge_gain_le_0 =
      sum(setting_summary$mean_ridge_gain_vs_intercept <= 0, na.rm = TRUE),
    min_mean_ridge_gain_vs_intercept =
      min(setting_summary$mean_ridge_gain_vs_intercept),
    median_mean_ridge_gain_vs_intercept =
      stats::median(setting_summary$mean_ridge_gain_vs_intercept),
    max_mean_ridge_gain_vs_intercept =
      max(setting_summary$mean_ridge_gain_vs_intercept),
    worst_setting_id =
      as.integer(worst$setting_id),
    worst_setting_a =
      as.numeric(worst$a),
    worst_setting_k =
      as.integer(worst$k),
    worst_setting_snr =
      as.numeric(worst$snr),
    worst_setting_mean_ridge_gain =
      as.numeric(worst$mean_ridge_gain_vs_intercept),
    worst_setting_ridge_gain_ci_lower =
      as.numeric(worst$ridge_gain_ci_lower),
    worst_setting_ridge_gain_ci_upper =
      as.numeric(worst$ridge_gain_ci_upper),
    min_mean_ridge_fraction_oracle_gain =
      min(setting_summary$mean_ridge_fraction_oracle_gain, na.rm = TRUE),
    median_mean_ridge_fraction_oracle_gain =
      stats::median(
        setting_summary$mean_ridge_fraction_oracle_gain,
        na.rm = TRUE
      ),
    min_prop_ridge_beats_intercept =
      min(setting_summary$prop_ridge_beats_intercept, na.rm = TRUE),
    analysis_note = paste(
      "Targeted final learnability diagnostic.",
      "Ridge is compared only with intercept/oracle.",
      "No Ridge-vs-LASSO or other between-model contrast is computed.",
      "CIs are descriptive."
    ),
    stringsAsFactors = FALSE
  ))
}

# Normalize batchtools job IDs returned as a table or vector.
# Inputs: x. Returns an integer vector, possibly empty.
ridge_intercept_job_ids <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(integer())
  }

  if (is.data.frame(x)) {
    if ("job.id" %in% names(x)) {
      return(as.integer(x$job.id))
    }
    if ("job.id" %in% colnames(x)) {
      return(as.integer(x[, "job.id"]))
    }
  }

  return(as.integer(x))
}

# Combine worker rows while retaining short error records.
# Inputs: result_list. Returns one data frame with absent fields filled by NA.
ridge_intercept_bind_results <- function(result_list) {
  good <- Filter(Negate(is.null), result_list)

  if (length(good) == 0L) {
    return(data.frame())
  }

  rows <- lapply(good, function(x) {
    if (is.data.frame(x)) {
      return(x)
    }
    return(as.data.frame(x, stringsAsFactors = FALSE))
  })

  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))

  rows <- lapply(rows, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) {
      d[[nm]] <- NA
    }
    return(d[, all_names, drop = FALSE])
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}

# Prepare a socket-worker script that loads the actual project definitions.
# Inputs: destination, project_root and source paths. Returns the written path.
# The generated script fixes BLAS threads to one and resolves sources absolutely.
ridge_intercept_write_worker_bootstrap <- function(
    destination,
    project_root = ".",
    pilot_core_file = "R/pilot_core.R",
    diagnostic_file = "R/ridge_intercept_check.R"
) {
  project_root <- normalizePath(
    project_root,
    winslash = "/",
    mustWork = TRUE
  )

  pilot_abs <- normalizePath(
    file.path(project_root, pilot_core_file),
    winslash = "/",
    mustWork = TRUE
  )

  diagnostic_abs <- normalizePath(
    file.path(project_root, diagnostic_file),
    winslash = "/",
    mustWork = TRUE
  )

  destination <- path.expand(destination)

  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  q <- function(x) return(encodeString(x, quote = '"'))

  lines <- c(
    "# Auto-generated Ridge/intercept batchtools worker bootstrap.",
    "# Source the real project files via absolute paths.",
    'Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",',
    '           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")',
    paste0("source(", q(pilot_abs), ")"),
    paste0("source(", q(diagnostic_abs), ")")
  )

  writeLines(
    lines,
    con = destination,
    useBytes = TRUE
  )

  if (!file.exists(destination)) {
    stop(
      "Failed to create batchtools worker bootstrap: ",
      destination,
      call. = FALSE
    )
  }

  return(normalizePath(
    destination,
    winslash = "/",
    mustWork = TRUE
  ))
}

# Record the setup and file hashes needed to validate a resumed registry.
# Inputs: run options, design tables and source_files. Returns a signature list.
ridge_intercept_signature <- function(
    R,
    n_test,
    n_cores,
    seed_base,
    full_settings,
    selected_settings,
    source_files
) {

  source_files <- unique(c(source_files, reproducibility_script_file))
  source_files <- normalizePath(
    source_files,
    winslash = "/",
    mustWork = TRUE
  )

  md5 <- as.character(tools::md5sum(source_files))
  names(md5) <- basename(source_files)

  tf_full <- tempfile(fileext = ".rds")
  tf_sel <- tempfile(fileext = ".rds")

  on.exit(
    unlink(c(tf_full, tf_sel)),
    add = TRUE
  )

  saveRDS(
    full_settings,
    tf_full,
    version = 2
  )
  saveRDS(
    selected_settings,
    tf_sel,
    version = 2
  )

  return(list(
    version = ridge_intercept_check_version,
    model = "Ridge_alpha_0",
    alpha = 0,
    R = as.integer(R),
    n_test = as.integer(n_test),
    seed_base = as.integer(seed_base),
    n_full_settings = nrow(full_settings),
    n_selected_settings = nrow(selected_settings),
    selected_setting_id =
      as.integer(selected_settings$setting_id),
    full_settings_md5 =
      unname(tools::md5sum(tf_full)),
    selected_settings_md5 =
      unname(tools::md5sum(tf_sel)),
    source_md5 = md5,
    cluster_function = "socket",
    n_cores = as.integer(n_cores)
  ))
}

# Validate task identity, model labels and successful worker rows.
# Inputs: ans, expected_task_id and model_spec. Returns invisible TRUE or stops.
ridge_intercept_validate_success_result <- function(
    ans,
    expected_task_id
) {
  if (!is.data.frame(ans)) {
    stop("Worker result is not a data.frame.", call. = FALSE)
  }

  if (nrow(ans) != 1L) {
    msg <- if ("error_message" %in% names(ans)) {
      paste(
        unique(na.omit(ans$error_message)),
        collapse = " | "
      )
    } else {
      "no error message"
    }

    stop(
      "Expected one Ridge result row for task ",
      expected_task_id,
      ", but got ",
      nrow(ans),
      ". ",
      msg,
      call. = FALSE
    )
  }

  if (!isTRUE(ans$ok[[1L]])) {
    stop(
      "Ridge worker result has ok != TRUE for task ",
      expected_task_id,
      ": ",
      if ("error_message" %in% names(ans)) ans$error_message[[1L]] else "unknown",
      call. = FALSE
    )
  }

  if (as.integer(ans$task_id[[1L]]) != as.integer(expected_task_id)) {
    stop("Worker task_id mismatch.", call. = FALSE)
  }

  if (!identical(as.character(ans$model[[1L]]), "Ridge_alpha_0")) {
    stop("Worker model label is not Ridge_alpha_0.", call. = FALSE)
  }

  if (!isTRUE(abs(as.numeric(ans$alpha[[1L]])) < 1e-15)) {
    stop("Worker alpha is not 0.", call. = FALSE)
  }

  return(invisible(TRUE))
}

# Execute one real task in a temporary socket registry before a large run.
# Inputs: design tables, task_row and worker options. Returns the checked result.
# Uses a smaller test sample and removes its temporary registry afterwards.
ridge_intercept_batchtools_smoke_test <- function(
    full_settings,
    selected_settings,
    design,
    task_row,
    analysis_n_test = 500L,
    project_root = ".",
    n_cores = 1L
) {
  ridge_intercept_assert_master_infrastructure()

  if (!requireNamespace("batchtools", quietly = TRUE)) {
    stop(
      "Package 'batchtools' is required for the socket smoke test.",
      call. = FALSE
    )
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "Package 'glmnet' is required for the socket smoke test.",
      call. = FALSE
    )
  }

  td <- tempfile(pattern = "ridge_intercept_bt_smoke_")

  dir.create(
    td,
    recursive = TRUE,
    showWarnings = FALSE
  )

  on.exit(
    unlink(td, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  bootstrap <- ridge_intercept_write_worker_bootstrap(
    destination = file.path(td, "worker_bootstrap.R"),
    project_root = project_root
  )

  reg <- batchtools::makeRegistry(
    file.dir = file.path(td, "registry"),
    packages = "glmnet",
    source = bootstrap,
    seed = 1L,
    make.default = FALSE
  )

  batchtools::batchMap(
    fun = ridge_intercept_bt_worker,
    task_row = list(
      as.data.frame(
        task_row,
        stringsAsFactors = FALSE
      )
    ),
    more.args = list(
      full_settings = full_settings,
      selected_settings = selected_settings,
      design = design,
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

  try(
    batchtools::syncRegistry(reg = reg),
    silent = TRUE
  )

  err <- ridge_intercept_job_ids(
    batchtools::findErrors(reg = reg)
  )
  not_done <- ridge_intercept_job_ids(
    batchtools::findNotDone(reg = reg)
  )

  if (length(err) > 0L || length(not_done) > 0L) {
    stop(
      "Batchtools/socket smoke test failed before the 600-job run. errors=",
      length(err),
      ", not_done=",
      length(not_done),
      ". The full run was NOT started.",
      call. = FALSE
    )
  }

  ans_list <- batchtools::reduceResultsList(
    reg = reg,
    missing.val = NULL
  )

  if (length(ans_list) != 1L || is.null(ans_list[[1L]])) {
    stop(
      "Batchtools/socket smoke test did not return exactly one task result.",
      call. = FALSE
    )
  }

  expected_task_id <- as.integer(
    as.data.frame(
      task_row,
      stringsAsFactors = FALSE
    )$task_id[1L]
  )

  ridge_intercept_validate_success_result(
    ans = ans_list[[1L]],
    expected_task_id = expected_task_id
  )

  return(invisible(ans_list[[1L]]))
}

# Run or resume the 12-setting Ridge check with 50 repetitions and 5000 test rows.
# Inputs: worker/seed options, output_dir, registry_dir and restart/save flags.
# Returns the registry, tasks, checked results and summaries invisibly.
# save_outputs writes tables, the summary object and session information.
# Failure rows are saved before stopping. This check computes no Ridge-LASSO contrast.
run_ridge_intercept_check <- function(
    R = 50L,
    n_test = 5000L,
    n_cores = 8L,
    seed_base = 1700000000L,
    output_dir = "outputs/ridge_intercept_check",
    registry_dir = file.path(output_dir, "registry_R50_corrected"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = TRUE,
    wait = TRUE,
    save_outputs = TRUE
) {

  ridge_intercept_assert_scientific_core()
  ridge_intercept_assert_master_infrastructure()

  if (!requireNamespace("batchtools", quietly = TRUE)) {
    stop("Package 'batchtools' is required.", call. = FALSE)
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required.", call. = FALSE)
  }

  R <- as.integer(R)
  n_test <- as.integer(n_test)
  n_cores <- as.integer(n_cores)
  seed_base <- as.integer(seed_base)

  if (R != 50L) {
    stop("This reviewed follow-up must use R=50.", call. = FALSE)
  }

  if (n_test != 5000L) {
    stop(
      "This reviewed follow-up must use analysis n_test=5000.",
      call. = FALSE
    )
  }

  if (n_cores < 1L) {
    stop("n_cores must be >= 1.", call. = FALSE)
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  design <- ridge_intercept_prepare_design()

  full_settings <- ridge_intercept_full_settings(
    design,
    expected_n_settings = 192L
  )

  selected_settings <- ridge_intercept_select_settings(
    full_settings
  )

  tasks <- ridge_intercept_make_tasks(
    selected_settings = selected_settings,
    R = R,
    seed_base = seed_base
  )

  if (nrow(tasks) != 600L) {
    stop(
      "Expected exactly 600 batch tasks (12 settings x 50 repetitions).",
      call. = FALSE
    )
  }

  signature_source_files <- c(
    "R/pilot_core.R",
    "R/precision_pilot.R",
    "R/ridge_intercept_check.R"
  )

  missing_source <- signature_source_files[
    !file.exists(signature_source_files)
  ]

  if (length(missing_source) > 0L) {
    stop(
      "Missing source file(s): ",
      paste(missing_source, collapse = ", "),
      call. = FALSE
    )
  }

  project_root_abs <- normalizePath(
    ".",
    winslash = "/",
    mustWork = TRUE
  )

  worker_bootstrap <- ridge_intercept_write_worker_bootstrap(
    destination = file.path(
      output_dir,
      "ridge_intercept_worker_bootstrap.R"
    ),
    project_root = project_root_abs
  )

  signature <- ridge_intercept_signature(
    R = R,
    n_test = n_test,
    n_cores = n_cores,
    seed_base = seed_base,
    full_settings = full_settings,
    selected_settings = selected_settings,
    source_files = signature_source_files
  )

  registry_exists <-
    dir.exists(registry_dir) &&
      file.exists(file.path(registry_dir, "registry.rds"))

  signature_file <- file.path(
    registry_dir,
    "ridge_intercept_signature.rds"
  )

  if (registry_exists && overwrite_registry) {
    unlink(
      registry_dir,
      recursive = TRUE,
      force = TRUE
    )
    registry_exists <- FALSE
  }

  if (registry_exists) {
    if (!resume) {
      stop(
        "Registry already exists and resume=FALSE: ",
        registry_dir,
        call. = FALSE
      )
    }

    if (!file.exists(signature_file)) {
      stop(
        "Existing registry has no signature file, refusing unsafe resume.",
        call. = FALSE
      )
    }

    old_signature <- readRDS(signature_file)

    if (!identical(old_signature, signature)) {
      stop(
        "Existing registry signature does not match this run. ",
        "Use a new registry path, do not overwrite completed scientific results.",
        call. = FALSE
      )
    }

    reg <- batchtools::loadRegistry(
      file.dir = registry_dir,
      writeable = TRUE,
      make.default = FALSE
    )

    try(
      batchtools::syncRegistry(reg = reg),
      silent = TRUE
    )
  } else {
    dir.create(
      dirname(registry_dir),
      recursive = TRUE,
      showWarnings = FALSE
    )

    reg <- batchtools::makeRegistry(
      file.dir = registry_dir,
      packages = "glmnet",
      source = worker_bootstrap,
      seed = 1L,
      make.default = FALSE
    )

    task_list <- split(
      tasks,
      seq_len(nrow(tasks))
    )

    batchtools::batchMap(
      fun = ridge_intercept_bt_worker,
      task_row = task_list,
      more.args = list(
        full_settings = full_settings,
        selected_settings = selected_settings,
        design = design,
        analysis_n_test = n_test
      ),
      reg = reg
    )

    dir.create(
      registry_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    saveRDS(
      signature,
      signature_file,
      version = 2
    )

    saveRDS(
      full_settings,
      file.path(registry_dir, "full_settings.rds"),
      version = 2
    )

    saveRDS(
      selected_settings,
      file.path(registry_dir, "selected_settings.rds"),
      version = 2
    )

    saveRDS(
      tasks,
      file.path(registry_dir, "tasks.rds"),
      version = 2
    )
  }

  reg <- set_batchtools_cluster_function(
    reg = reg,
    n_cores = n_cores,
    cluster_function = "socket"
  )

  on_system <- ridge_intercept_job_ids(
    batchtools::findOnSystem(reg = reg)
  )

  if (length(on_system) > 0L) {
    stop(
      "Registry still reports ",
      length(on_system),
      " job(s) as 'on system'. Do not start a second copy. ",
      "After a confirmed server restart, stale jobs can be reset before resuming.",
      call. = FALSE
    )
  }

  err_ids <- ridge_intercept_job_ids(
    batchtools::findErrors(reg = reg)
  )

  if (length(err_ids) > 0L && isTRUE(resubmit_errors)) {
    batchtools::resetJobs(
      ids = err_ids,
      reg = reg
    )
  }

  todo <- ridge_intercept_job_ids(
    batchtools::findNotDone(reg = reg)
  )

  if (length(todo) > 0L) {
    message(
      "Submitting ",
      length(todo),
      " unfinished Ridge/intercept jobs out of ",
      nrow(tasks),
      "."
    )

    batchtools::submitJobs(
      ids = todo,
      reg = reg
    )

    if (isTRUE(wait)) {
      batchtools::waitForJobs(
        ids = todo,
        reg = reg
      )
    }
  }

  if (!isTRUE(wait)) {
    return(
      invisible(
        list(
          registry = reg,
          full_settings = full_settings,
          selected_settings = selected_settings,
          tasks = tasks,
          signature = signature
        )
      )
    )
  }

  try(
    batchtools::syncRegistry(reg = reg),
    silent = TRUE
  )

  err_after <- ridge_intercept_job_ids(
    batchtools::findErrors(reg = reg)
  )

  if (length(err_after) > 0L && isTRUE(resubmit_errors)) {
    message(
      "Retrying ",
      length(err_after),
      " batchtools error job(s) once."
    )

    batchtools::resetJobs(
      ids = err_after,
      reg = reg
    )

    batchtools::submitJobs(
      ids = err_after,
      reg = reg
    )

    batchtools::waitForJobs(
      ids = err_after,
      reg = reg
    )

    try(
      batchtools::syncRegistry(reg = reg),
      silent = TRUE
    )
  }

  remaining <- ridge_intercept_job_ids(
    batchtools::findNotDone(reg = reg)
  )

  if (length(remaining) > 0L) {
    stop(
      length(remaining),
      " jobs are still not done after wait/retry. ",
      "Inspect the registry before resuming.",
      call. = FALSE
    )
  }

  result_list <- batchtools::reduceResultsList(
    reg = reg,
    missing.val = NULL
  )

  if (length(result_list) != nrow(tasks)) {
    stop(
      "Collected ",
      length(result_list),
      " task result objects but expected ",
      nrow(tasks),
      ".",
      call. = FALSE
    )
  }

  validation_errors <- character()

  for (i in seq_along(result_list)) {
    ans <- result_list[[i]]

    if (is.null(ans)) {
      validation_errors <- c(
        validation_errors,
        paste0("task ", tasks$task_id[i], ": NULL result")
      )
      next
    }

    validation_error <- tryCatch({
      ridge_intercept_validate_success_result(
        ans = ans,
        expected_task_id = tasks$task_id[i]
      )
      NULL

    }, error = function(e) {
      return(paste0("task ", tasks$task_id[i], ": ", conditionMessage(e)))
    })
    if (!is.null(validation_error)) {
      validation_errors <- c(validation_errors, validation_error)
    }
  }

  results <- ridge_intercept_bind_results(
    result_list
  )

  failure_records <- results[is.na(results$ok) | !results$ok, , drop = FALSE]
  if (nrow(failure_records) > 0L && isTRUE(save_outputs)) {
    utils::write.csv(
      failure_records, file.path(output_dir, "ridge_intercept_failures.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      results, file.path(output_dir, "ridge_intercept_task_results_partial.csv"),
      row.names = FALSE
    )
  }

  if (length(validation_errors) > 0L) {
    utils::write.csv(
      results,
      file.path(
        output_dir,
        "ridge_intercept_task_results_partial.csv"
      ),
      row.names = FALSE
    )

    writeLines(
      validation_errors,
      con = file.path(
        output_dir,
        "ridge_intercept_validation_errors.txt"
      )
    )

    stop(
      length(validation_errors),
      " task result(s) failed structural validation. ",
      "Partial outputs were saved, do not interpret them.",
      call. = FALSE
    )
  }

  if (nrow(results) != nrow(tasks)) {
    stop(
      "Collected ",
      nrow(results),
      " result rows but expected ",
      nrow(tasks),
      ".",
      call. = FALSE
    )
  }

  results <- results[
    order(results$task_id),
    ,
    drop = FALSE
  ]
  rownames(results) <- NULL

  if (!identical(
    as.integer(results$task_id),
    seq_len(nrow(tasks))
  )) {
    stop(
      "Task-result IDs are incomplete or out of alignment.",
      call. = FALSE
    )
  }

  failed <- results[
    is.na(results$ok) | !results$ok,
    ,
    drop = FALSE
  ]

  successful <- results[
    !is.na(results$ok) & results$ok,
    ,
    drop = FALSE
  ]

  if (nrow(successful) > 0L) {
    setting_summary <- ridge_intercept_setting_summary(
      successful
    )

    by_snr <- ridge_intercept_by_snr(
      setting_summary
    )

    overview <- ridge_intercept_overview(
      results = results,
      setting_summary = setting_summary,
      expected_jobs = nrow(tasks)
    )

    lowest_gain <- setting_summary[
      order(setting_summary$mean_ridge_gain_vs_intercept),
      ,
      drop = FALSE
    ]
  } else {
    setting_summary <- data.frame()
    by_snr <- data.frame()
    overview <- data.frame()
    lowest_gain <- data.frame()
  }

  if (isTRUE(save_outputs)) {
    utils::write.csv(
      results,
      file.path(
        output_dir,
        "ridge_intercept_task_results.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      selected_settings,
      file.path(
        output_dir,
        "ridge_intercept_selected_settings.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      setting_summary,
      file.path(
        output_dir,
        "ridge_intercept_setting_summary.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      by_snr,
      file.path(
        output_dir,
        "ridge_intercept_by_snr.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      overview,
      file.path(
        output_dir,
        "ridge_intercept_overview.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      lowest_gain,
      file.path(
        output_dir,
        "ridge_intercept_lowest_gain_settings.csv"
      ),
      row.names = FALSE
    )

    utils::write.csv(
      failed,
      file.path(
        output_dir,
        "ridge_intercept_failures.csv"
      ),
      row.names = FALSE
    )

    saveRDS(
      list(
        version = ridge_intercept_check_version,
        signature = signature,
        full_settings = full_settings,
        selected_settings = selected_settings,
        tasks = tasks,
        results = results,
        setting_summary = setting_summary,
        by_snr = by_snr,
        overview = overview,
        registry_status = batchtools::getStatus(reg = reg)
      ),
      file.path(
        output_dir,
        "ridge_intercept_summary.rds"
      ),
      version = 2
    )

    capture.output(
      sessionInfo(),
      file = file.path(
        output_dir,
        "ridge_intercept_sessionInfo.txt"
      )
    )
  }

  if (nrow(failed) > 0L) {
    stop(
      "The batchtools run finished, but ",
      nrow(failed),
      " task(s) returned ok=FALSE. ",
      "Outputs were saved, inspect ridge_intercept_failures.csv before interpretation.",
      call. = FALSE
    )
  }

  return(invisible(
    list(
      registry = reg,
      full_settings = full_settings,
      selected_settings = selected_settings,
      tasks = tasks,
      results = results,
      setting_summary = setting_summary,
      by_snr = by_snr,
      overview = overview,
      signature = signature
    )
  ))
}
