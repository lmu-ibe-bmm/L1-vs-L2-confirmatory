###############################################################
# BLINDED PRECISION PILOT FOR SELECTING THE REPETITION COUNT
###############################################################
#
# PURPOSE
# This script estimates the number of Monte Carlo repetitions required for
# the primary relative prediction mean squared error estimand comparing Ridge
# and the least absolute shrinkage and selection operator, abbreviated below as
# LASSO. Ridge is deliberately not fitted in this pilot, so the confirmatory
# comparison is not examined during planning. LASSO and Random Forest provide
# information about variability and paired dependence instead.
#
# CONTRIBUTION TO THE OVERALL STUDY
# The test-set-size pilot controls uncertainty caused by evaluating predictions
# on a finite test sample. This precision pilot addresses a different source of
# uncertainty, namely Monte Carlo variation across independently generated
# study repetitions. Its only decision is the number of repetitions needed for
# sufficient precision. It does not provide evidence for or against any
# Ridge-to-LASSO hypothesis.
#
# STATISTICAL IDEA
# For theta = 1 - E(MSE_C) / E(MSE_LASSO), where MSE denotes prediction mean
# squared error, planning is based on the standard deviation of the paired
# influence values. If tau denotes that standard deviation and R denotes the
# number of repetitions, the Monte Carlo standard error of the estimated theta
# is approximately tau divided by the square root of R.
#
# The Random Forest mean is never treated as a Ridge mean. Only centered Random
# Forest deviations are used. Their correlation with LASSO is preserved, while
# the comparator standard deviation is set to 1.0, 1.1, 1.2, or 1.3 times the
# LASSO standard deviation. These scenarios make the final repetition choice
# less dependent on the raw Random Forest scale.
#
# HOW THE PILOT OUTPUTS ARE USED
# Each task calculates paired LASSO and Random Forest prediction mean squared
# errors from the same simulated data. These paired values are required
# temporarily because their correlation affects the precision calculation.
# Within each setting, the code estimates the two standard deviations and their
# correlation. It then constructs the prespecified variance-only Random Forest
# surrogate for every standard-deviation multiplier and local-effect
# combination. The required number of repetitions is calculated from the
# influence-value standard deviation for each scenario. Finally, the code
# takes the largest requirement over scenarios and settings and rounds it
# upward to a block of one hundred repetitions. The public outputs contain the
# resulting standard deviations, coefficients of variation, dependence
# measures, and repetition requirements, not the individual prediction errors.
# This is a conservative planning rule, not a test of a scientific hypothesis.
#
# REPORT FIGURE
# In addition to the planning tables, this script creates a figure showing,
# for every setting,
#        SD(PMSE_RF) / SD(PMSE_LASSO),
# where SD means standard deviation and PMSE means prediction mean squared
# error.
# These descriptive raw ratios are not planning assumptions. Planning continues
# to use the separate sensitivity scenarios listed above.
#
# Each figure is tied to its input through a source tag. The helper
# resolve_precision_figure_source() accepts the latest timestamped setting
# summary or an explicit comma-separated values file path.
#
# METHODOLOGICAL BASIS
# The separation of aims, data generation, estimands, methods, and performance
# measures follows Morris et al. (2019) and the planning and reporting
# framework of Siepe et al. (2024). The ratio-of-means influence value uses the
# delta method described by Deng et al. (2018). Rearranging the Monte Carlo
# standard-error relationship gives R = ceiling((tau / target)^2).
#
# These references support the statistical structure. They do not prescribe
# this study's numerical target, minimum repetition count, Random Forest
# surrogate, sensitivity grid, or rounding rule. Those are prespecified study
# choices stored explicitly in the design and output tables.
#
# WHEN AND HOW TO RUN THIS SCRIPT
# Do not start it by sourcing this file interactively unless custom arguments
# are required. From the top-level package directory, first run
# run/00_install_packages.R, then tests/run_all_tests.R, and finally
# run/02_run_precision_pilot.R. The launcher writes a summary R data file,
# three readable summary tables, and timestamped Portable Network Graphics and
# Portable Document Format figures to outputs/precision_pilot/. Individual
# prediction errors remain only in the local batchtools registry while that
# registry is needed to resume or summarize an interrupted run. The registry is
# execution state. To plot a selected summary without rerunning the pilot,
# use run/04_export_precision_figure.R.
###############################################################

# PURPOSE
# Locate another R file stored in the same directory as this script.
#
# WHY THIS FUNCTION EXISTS
# A user may start a launcher from the top-level package directory, source this
# file from an interactive R session, or call it with Rscript from another
# directory.
# Relative paths behave differently in those situations. This helper searches
# the location of the running/sourced script before falling back to the current
# working directory, which makes the subsequent source("pilot_core.R") call
# independent of how the project was started.
#
# INPUT
# filename: File name to locate, for example "pilot_core.R".
#
# RETURN VALUE
# A normalized absolute path to the first existing candidate. The function
# stops with a clear error if the file cannot be found.
find_sibling_file <- function(filename) {
  # Candidate 1: the directory of a file passed to Rscript via --file=.
  candidates <- character(0)

  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", cmd, fixed = TRUE)
  if (length(hit) > 0L) {
    this_file <- sub("--file=", "", cmd[hit[1]], fixed = TRUE)
    this_file <- normalizePath(this_file, winslash = "/", mustWork = FALSE)
    candidates <- c(candidates, file.path(dirname(this_file), filename))
  }

  # Candidate 2: the directory recorded by source() in the active call stack.
  frames <- sys.frames()
  ofiles <- vapply(
    frames,
    function(x) return(if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_),
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) {
    this_file <- normalizePath(
      ofiles[length(ofiles)], winslash = "/", mustWork = FALSE
    )
    candidates <- c(candidates, file.path(dirname(this_file), filename))
  }

  # Candidate 3: the caller's current working directory.
  candidates <- unique(c(candidates, filename))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) stop("File not found: ", filename)

  # Returning an absolute path prevents later working-directory changes from
  # changing which file is sourced.
  return(normalizePath(existing[1], winslash = "/", mustWork = TRUE))
}

# Load the shared simulation, model-fitting, and execution functions. The core
# also fixes all basic linear algebra subprogram and model-internal thread counts
# before any computation,
# so those environment settings are not duplicated in this module.
script_file <- find_sibling_file("pilot_core.R")
source(script_file)


###############################################################
# A. SELECTING A PRECISION RESULT FOR THE FIGURE
###############################################################

# PURPOSE
# Select a precision setting-summary table and its source-specific figure name.
#
# WHY THIS FUNCTION EXISTS
# Explicit paths connect figures to one execution. The latest mode selects
# timestamped files by their names rather than filesystem modification times.
#
# INPUTS
# selection: Either latest or one existing CSV path.
# results_dir: Directory searched when selection is latest.
# base_figure_stem: Figure basename before the source tag is appended.
#
# RETURN VALUE
# A list with mode, input_file, source_tag and figure_file_stem. Missing inputs
# stop execution before a figure is written.
resolve_precision_figure_source <- function(
    selection = "latest",
    results_dir = file.path("outputs", "precision_pilot"),
    base_figure_stem = "Supp_Fig_RF_LASSO_SD_ratio_boxplot_pretty") {
  if (length(selection) != 1L || is.na(selection) || !nzchar(trimws(selection))) {
    stop("selection must be 'latest' or one existing CSV path.")
  }
  selection <- trimws(selection)
  mode <- if (tolower(selection) == "latest") "latest" else "file"
  if (mode == "latest") {
    pattern <- paste0(
      "^sd_precision_setting_summary_",
      "[0-9]{4}-[0-9]{2}-[0-9]{2}_",
      "[0-9]{2}-[0-9]{2}-[0-9]{2}\\.csv$"
    )
    candidates <- list.files(results_dir, pattern = pattern, full.names = TRUE)
    if (length(candidates) == 0L) stop("No precision setting summary in: ", results_dir)
    input_file <- tail(sort(candidates), 1L)
  } else {
    input_file <- path.expand(selection)
  }
  if (!file.exists(input_file)) stop("Precision setting summary not found: ", input_file)
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  tag <- sub("^sd_precision_setting_summary_", "",
             tools::file_path_sans_ext(basename(input_file)))
  tag <- gsub("[^A-Za-z0-9_-]+", "_", tag)
  return(list(
    mode = mode, input_file = input_file, source_tag = tag,
    figure_file_stem = paste0(base_figure_stem, "_", tag)
  ))
}


###############################################################
# B. REPORT FIGURE: SETTING-SPECIFIC PREDICTION-ERROR VARIABILITY
###############################################################

# PURPOSE
# Convert a setting-level precision summary into the data and annotations used
# by the Random-Forest/LASSO variability figure.
#
# WHY THIS FUNCTION EXISTS
# Figure preparation is kept separate from drawing and saving the plot. This
# makes the numerical values testable and ensures that the displayed counts
# are calculated from the supplied simulation output rather than typed into
# the plotting code.
#
# INPUT
# setting_summary: One row per simulation setting. It must contain setting_id,
#   sd_mse_lasso, and sd_mse_rf_raw.
#
# RETURN VALUE
# A list with:
#   data    - the three plotting columns, including the Random Forest to LASSO
#             standard-deviation ratio,
#   summary - the median ratio and counts above the reference thresholds.
prepare_rf_lasso_sd_ratio <- function(setting_summary) {
  # First confirm that the table contains the quantities needed for the ratio.
  required_cols <- c("setting_id", "sd_mse_lasso", "sd_mse_rf_raw")
  missing_cols <- setdiff(required_cols, names(setting_summary))

  if (length(missing_cols) > 0L) {
    stop(
      paste(
        "Columns required for the Random Forest to LASSO",
        "standard-deviation figure are missing:"
      ),
      " ",
      paste(missing_cols, collapse = ", "), "."
    )
  }

  # Keep a small, self-contained plotting data set and give the Random Forest
  # column the shorter field name used in the remainder of this function.
  plot_data <- setting_summary[, required_cols, drop = FALSE]
  names(plot_data)[names(plot_data) == "sd_mse_rf_raw"] <- "sd_mse_rf"

  if (any(!is.finite(plot_data$sd_mse_lasso)) ||
    any(!is.finite(plot_data$sd_mse_rf)) ||
    any(plot_data$sd_mse_lasso <= 0) ||
    any(plot_data$sd_mse_rf < 0)) {
    stop(
      "The setting-specific prediction mean squared error standard deviations ",
      "must be finite and ",
      "non-negative, with positive LASSO standard deviations."
    )
  }

  # Each ratio compares variability across repetitions within the same setting.
  plot_data$sd_ratio_rf_lasso <-
    plot_data$sd_mse_rf / plot_data$sd_mse_lasso

  # These values are printed in the figure. They are calculated from the data
  # rather than hard-coded, so a rerun automatically updates the annotations.
  figure_summary <- data.frame(
    n_settings = nrow(plot_data),
    median_ratio = stats::median(plot_data$sd_ratio_rf_lasso),
    n_above_1 = sum(plot_data$sd_ratio_rf_lasso > 1),
    n_above_1_3 = sum(plot_data$sd_ratio_rf_lasso > 1.3),
    n_above_2 = sum(plot_data$sd_ratio_rf_lasso > 2),
    stringsAsFactors = FALSE
  )

  return(list(data = plot_data, summary = figure_summary))
}


# PURPOSE
# Draw and save the report figure comparing setting-specific Random Forest and
# LASSO prediction mean squared error standard deviations.
#
# WHY THIS FUNCTION EXISTS
# The visual specification lives next to the calculation that supplies it, so
# rerunning the pilot or the figure-only launcher reproduces the same graph.
# The plot is descriptive. Its raw Random Forest to LASSO ratios are not
# substituted for the sensitivity assumptions used to choose the planned
# number of repetitions.
#
# INPUTS
# setting_summary: Setting-level table accepted by
#   prepare_rf_lasso_sd_ratio().
# output_dir: Directory in which the Portable Network Graphics and Portable
#   Document Format files are written.
# file_stem: Shared file name without an extension. The direct-call default is
#   timestamped to keep figures from separate executions distinguishable.
# width_in, height_in: Physical dimensions in inches.
# dpi: Dots-per-inch resolution of the Portable Network Graphics output.
#
# RETURN VALUE
# A list containing both file paths, the numerical annotation summary, and the
# ggplot object. Returning the plot object allows inspection without rereading
# an image file.
save_rf_lasso_sd_ratio_figure <- function(
    setting_summary,
    output_dir,
    file_stem = paste0(
      "Supp_Fig_RF_LASSO_SD_ratio_boxplot_pretty_",
      make_timestamp()
    ),
    width_in = 7.1,
    height_in = 5.3,
    dpi = 300) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for the report figure. ",
      "Install it and rerun the plotting step."
    )
  }

  # Step 1: calculate the ratio and all text annotations from the supplied data.
  prepared <- prepare_rf_lasso_sd_ratio(setting_summary)
  d <- prepared$data
  s <- prepared$summary

  # A one-level factor produces the single box shown on the x-axis.
  d$plot_group <- factor(
    "RF SD / LASSO SD",
    levels = "RF SD / LASSO SD"
  )

  # Assemble the upper-right annotation from calculated, not hard-coded, values.
  summary_label <- paste0(
    "Median = ",
    formatC(s$median_ratio, format = "f", digits = 2, decimal.mark = "."),
    "\n",
    "> 1.0 in ", s$n_above_1, "/", s$n_settings, " settings\n",
    "> 1.3 in ", s$n_above_1_3, "/", s$n_settings, " settings\n",
    "> 2.0 in ", s$n_above_2, "/", s$n_settings, " settings"
  )

  # Step 2: assemble the report figure from the calculated data and labels.
  # The layers reproduce a horizontal grey grid, blue reference lines, a black
  # boxplot, an orange median line, and four setting-level summary statistics
  # in the upper-right corner. Explicit moderate text sizes keep every label
  # inside the plotting panel across the graphics devices used on the server
  # and on ordinary desktop installations.
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = plot_group, y = sd_ratio_rf_lasso)
  ) +
    ggplot2::geom_hline(
      yintercept = seq(0, 3, by = 0.5),
      colour = "#B0B0B0", linewidth = 0.45, alpha = 0.38
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      colour = "#1F77B4", linewidth = 0.75, linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 1.3,
      colour = "#1F77B4", linewidth = 0.75, linetype = "dotted"
    ) +
    ggplot2::geom_boxplot(
      width = 0.45,
      fill = "white",
      colour = "black",
      linewidth = 0.75,
      outlier.shape = 21,
      outlier.fill = "white",
      outlier.colour = "black",
      outlier.stroke = 0.8,
      outlier.size = 2.5
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.775, xend = 1.225,
      y = s$median_ratio, yend = s$median_ratio,
      colour = "#FF7F0E",
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text", x = 1.02, y = 1.04,
      label = "equal variability", hjust = 0, vjust = 0,
      size = 3.5
    ) +
    ggplot2::annotate(
      "text", x = 1.02, y = 1.34,
      label = "30% higher RF variability", hjust = 0, vjust = 0,
      size = 3.5
    ) +
    ggplot2::annotate(
      "text", x = 1.15, y = 2.62,
      label = summary_label, hjust = 0, vjust = 1,
      size = 3.5, lineheight = 1.08
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 3.2),
      breaks = seq(0, 3, by = 0.5),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(
      x = NULL,
      y = "Ratio of setting-specific SD(PMSE)"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.75
      ),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.title.y = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )

  # Step 3: create the destination and save the same plot in a raster format
  # for convenient viewing and a vector format for the manuscript.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_file <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_stem, ".pdf"))

  ggplot2::ggsave(
    png_file, plot = p,
    width = width_in, height = height_in, units = "in", dpi = dpi,
    bg = "white"
  )
  # Select Cairo PDF when enabled, otherwise use the base R PDF device.
  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    pdf_file, plot = p,
    width = width_in, height = height_in, units = "in",
    device = pdf_device, bg = "white"
  )

  return(list(
    png = png_file,
    pdf = pdf_file,
    summary = s,
    plot = p
  ))
}

###############################################################
# C. VARIANCE-ONLY SURROGATE FOR PRECISION PLANNING
###############################################################
# PURPOSE
# Construct a controlled comparator vector from the observed Random Forest
# deviations for precision planning.
#
# WHY THIS FUNCTION EXISTS
# Random Forest is used only to learn a plausible paired dependence pattern. It
# is not a proxy estimate of Ridge's mean performance. Centering removes the raw
# Random Forest mean, scaling sets its standard deviation to a chosen multiple
# of the LASSO standard deviation, and shifting imposes a local assumed effect.
# These operations separate the correlation, variability, and mean assumptions.
#
# INPUTS
# mse_lasso: Paired LASSO prediction mean squared errors from one setting.
# mse_rf: Paired Random Forest prediction mean squared errors from the same
#   repetitions.
# sd_multiplier: Desired comparator standard deviation divided by the observed
#   LASSO standard deviation.
# theta_assumed: Assumed local value of 1 - mean(comparator)/mean(LASSO).
#
# RETURN VALUE
# A list containing the constructed comparator values, the resulting
# ratio-of-means precision summary, raw and surrogate moments, correlations,
# and the assumptions used to create the surrogate.
make_variance_only_rf_surrogate <- function(
    mse_lasso,
    mse_rf,
    sd_multiplier = 1,
    theta_assumed = 0) {

  # Step 1: retain only repetitions for which both prediction mean squared
  # errors are valid.
  ok <- is.finite(mse_lasso) & is.finite(mse_rf)
  L <- as.numeric(mse_lasso[ok])
  rf_mse <- as.numeric(mse_rf[ok])

  if (length(L) < 2L) stop("At least two valid pilot repetitions are required.")
  if (!is.finite(sd_multiplier) || sd_multiplier <= 0) {
    stop("sd_multiplier must be positive and finite.")
  }

  # Step 2: calculate the raw moments needed for centering and rescaling.
  mu_L <- mean(L)
  mu_F <- mean(rf_mse)
  sd_L <- sd(L)
  sd_F <- sd(rf_mse)

  if (!is.finite(mu_L) || mu_L <= 0 || !is.finite(sd_L) || sd_L < 0) {
    stop("The LASSO pilot moments are not valid.")
  }
  if (!is.finite(sd_F) || sd_F < 0) {
    stop("The Random Forest pilot standard deviation is not valid.")
  }

  # Step 3: translate the relative standard-deviation scenario into the desired
  # absolute standard deviation.
  target_sd <- sd_multiplier * sd_L

  # Step 4: use Random Forest deviations only for their correlation pattern.
  # Their standard deviation is deliberately replaced by the value defined in
  # the sensitivity scenario.
  if (sd_F > 0) {
    scale_factor <- target_sd / sd_F
  } else if (target_sd == 0) {
    scale_factor <- 0
  } else {
    stop(
      paste(
        "The Random Forest standard deviation is zero, so the requested",
        "positive surrogate standard deviation cannot be created."
      )
    )
  }

  # Step 5: center the Random Forest deviations and impose the assumed
  # comparator mean. Centering separates variability from the often very
  # different Random Forest mean.
  F_centered_scaled <- scale_factor * (rf_mse - mu_F)

  # theta = 1 - mu_C/mu_L  =>  mu_C = (1-theta)*mu_L.
  target_mean <- (1 - theta_assumed) * mu_L
  F_star <- target_mean + F_centered_scaled

  if (any(!is.finite(F_star))) {
    stop("The surrogate contains non-finite values.")
  }
  if (mean(F_star) <= 0) {
    stop("The surrogate mean must be positive.")
  }

  # Step 6: calculate the influence-value standard deviation tau. This is the
  # quantity used to derive the required number of Monte Carlo repetitions.
  ss <- ratio_of_means_summary(L, F_star)

  return(list(
    mse_surrogate = F_star,
    summary = ss,
    mu_lasso = mu_L,
    mu_rf_raw = mu_F,
    sd_lasso = sd_L,
    sd_rf_raw = sd_F,
    sd_rf_surrogate = sd(F_star),
    correlation_lasso_rf_raw = cor(L, rf_mse),
    correlation_lasso_surrogate = cor(L, F_star),
    sd_multiplier = sd_multiplier,
    scale_factor = scale_factor,
    theta_assumed = theta_assumed,
    target_mean_surrogate = target_mean,
    target_sd_surrogate = target_sd
  ))
}

###############################################################
# D. SETTING-LEVEL SUMMARY AND GLOBAL REPETITION SELECTION
###############################################################
# PURPOSE
# Summarize the completed precision pilot and select a conservative global
# repetition count for the planned study.
#
# WHY THIS FUNCTION EXISTS
# Precision requirements can differ across simulation settings. Therefore the
# function evaluates every setting under every prespecified combination of
# comparator standard-deviation multiplier and local effect assumption. The
# final planning value is the largest requirement across the full design and is
# rounded up to a convenient block of one hundred repetitions.
# Taking the maximum protects the requested precision in the most demanding
# examined scenario. Rounding upward is an operational choice that simplifies
# workload planning and cannot make the resulting Monte Carlo standard error
# larger than the value implied by the unrounded requirement. Neither operation
# is presented as a general recommendation from Morris et al. (2019).
#
# INPUTS
# run_results: One row per completed pilot repetition, including paired LASSO
#   and Random Forest prediction mean squared errors.
# settings: The setting table created by make_all_settings().
# design: Shared constants, including the Monte Carlo standard-error target and
#   minimum repetition count.
# sd_multiplier_grid: Relative comparator standard-deviation scenarios.
# theta_assumed_grid: Local effect values used to set the surrogate mean.
#
# RETURN VALUE
# A list with three tables. setting_summary records method-specific standard
# deviations and coefficients of variation, paired dependence, and the most
# demanding sensitivity scenario in each setting. It does not contain
# repetition-level prediction errors or method-specific mean prediction
# errors. sensitivity_summary records the global requirement separately for
# every comparator standard-deviation multiplier, which reproduces the
# sensitivity table in the Registered Report. global_summary records the final
# conservative maximum across all settings and sensitivity assumptions.
summarise_precision_pilot <- function(
    run_results,
    settings,
    design = make_design(),
    sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3),
    theta_assumed_grid = c(-0.02, 0, 0.02)) {

  # Step 1: confirm that all quantities needed for paired precision planning
  # are present and that no failed or duplicated task enters the summary.
  required_cols <- c(
    "task_id", "setting_id", "run_id", "ok", "mse_lasso", "mse_rf"
  )
  missing_cols <- setdiff(required_cols, names(run_results))
  if (length(missing_cols) > 0L) {
    stop(
      "The precision-pilot results are missing required columns: ",
      paste(missing_cols, collapse = ", "), "."
    )
  }

  if (any(!run_results$ok)) {
    stop("The precision pilot contains failed tasks and cannot be summarized.")
  }
  if (anyDuplicated(run_results$task_id) > 0L) {
    stop("The precision-pilot results contain duplicate task_id values.")
  }

  # Step 2: process settings separately because both the paired correlation
  # and the prediction mean squared error scales are setting-specific.
  ids <- sort(unique(as.integer(run_results$setting_id)))
  good <- run_results
  rows <- vector("list", length(ids))
  sensitivity_details <- vector("list", length(ids))

  for (ii in seq_along(ids)) {
    id <- ids[ii]
    d <- good[good$setting_id == id, , drop = FALSE]

    # A standard deviation cannot be estimated from fewer than two repetitions.
    # Preserve the setting in the output, but mark its precision fields missing.
    if (nrow(d) < 2L) {
      sensitivity_details[[ii]] <- tibble()
      rows[[ii]] <- tibble(
        setting_id = id,
        n_success = nrow(d),
        sd_mse_lasso = safe_sd(d$mse_lasso),
        cv_mse_lasso = NA_real_,
        sd_mse_rf_raw = safe_sd(d$mse_rf),
        cv_mse_rf_raw = NA_real_,
        corr_lasso_rf_raw = NA_real_,
        raw_tau_rf_vs_lasso = NA_real_,
        raw_R_required = NA_integer_,
        surrogate_scale_factor = NA_real_,
        surrogate_sd_multiplier_at_max = NA_real_,
        surrogate_sd_mse_rf = NA_real_,
        surrogate_corr_lasso_rf = NA_real_,
        surrogate_tau_max = NA_real_,
        surrogate_theta_assumed_at_max = NA_real_,
        R_required_sensitivity = NA_integer_
      )
      next
    }

    # The raw Random Forest comparison is retained for descriptive reporting
    # only.
    raw <- ratio_of_means_summary(d$mse_lasso, d$mse_rf)

    # Combine each standard-deviation multiplier with every local theta
    # assumption.
    sensitivity_rows <- bind_rows(lapply(sd_multiplier_grid, function(mult) {
      return(bind_rows(lapply(theta_assumed_grid, function(th) {
        z <- make_variance_only_rf_surrogate(
          mse_lasso = d$mse_lasso,
          mse_rf = d$mse_rf,
          sd_multiplier = mult,
          theta_assumed = th
        )

        return(tibble(
          sd_multiplier = mult,
          theta_assumed = th,
          tau = z$summary$tau_influence,
          mcse_at_Rpilot = z$summary$mcse,
          scale_factor = z$scale_factor,
          sd_rf_surrogate = z$sd_rf_surrogate,
          corr_surrogate = z$correlation_lasso_surrogate
        ))
      })))
    }))

    # Store every sensitivity scenario before selecting a setting-specific
    # maximum. These rows are needed to report one global result for each
    # comparator standard-deviation multiplier rather than only the final
    # multiplier selected for planning.
    sensitivity_rows$setting_id <- id
    sensitivity_rows$R_required <- vapply(
      sensitivity_rows$tau,
      required_R_from_tau,
      integer(1),
      mcse_target = design$mcse_target,
      R_min = design$R_min
    )
    sensitivity_details[[ii]] <- sensitivity_rows

    # The largest tau gives the largest repetition count for a fixed Monte Carlo
    # standard-error target. Retain the corresponding scenario so the
    # conservative choice remains auditable.
    idx_max <- which.max(sensitivity_rows$tau)
    max_row <- sensitivity_rows[idx_max, , drop = FALSE]

    rows[[ii]] <- tibble(
      setting_id = id,
      n_success = nrow(d),
      sd_mse_lasso = sd(d$mse_lasso),
      cv_mse_lasso = sd(d$mse_lasso) / mean(d$mse_lasso),
      sd_mse_rf_raw = sd(d$mse_rf),
      cv_mse_rf_raw = sd(d$mse_rf) / mean(d$mse_rf),
      corr_lasso_rf_raw = cor(d$mse_lasso, d$mse_rf),
      raw_tau_rf_vs_lasso = raw$tau_influence,
      raw_R_required = required_R_from_tau(
        raw$tau_influence,
        design$mcse_target,
        design$R_min
      ),
      surrogate_scale_factor = max_row$scale_factor,
      surrogate_sd_multiplier_at_max = max_row$sd_multiplier,
      surrogate_sd_mse_rf = max_row$sd_rf_surrogate,
      surrogate_corr_lasso_rf = max_row$corr_surrogate,
      surrogate_tau_max = max_row$tau,
      surrogate_theta_assumed_at_max = max_row$theta_assumed,
      R_required_sensitivity = required_R_from_tau(
        max_row$tau,
        design$mcse_target,
        design$R_min
      )
    )
  }

  # Step 3: attach the original factor values to the calculated summaries.
  setting_summary <- bind_rows(rows) %>%
    left_join(settings, by = "setting_id")

  # Step 4: summarize the most demanding setting and local-effect assumption
  # separately for every standard-deviation multiplier. This produces the four
  # rows reported for multipliers 1.0, 1.1, 1.2, and 1.3.
  sensitivity_by_setting <- bind_rows(sensitivity_details)
  sensitivity_summary <- bind_rows(lapply(
    sort(unique(as.numeric(sd_multiplier_grid))),
    function(mult) {
      d <- sensitivity_by_setting[
        sensitivity_by_setting$sd_multiplier == mult &
          is.finite(sensitivity_by_setting$tau) &
          is.finite(sensitivity_by_setting$R_required),
        ,
        drop = FALSE
      ]

      if (nrow(d) == 0L) {
        return(tibble(
          sd_multiplier = mult,
          max_tau = NA_real_,
          setting_id_at_max = NA_integer_,
          theta_assumed_at_max = NA_real_,
          R_global_unrounded = NA_integer_,
          R_global_rounded_100 = NA_integer_
        ))
      }

      index_max <- which.max(d$tau)
      global_required <- max(d$R_required)
      return(tibble(
        sd_multiplier = mult,
        max_tau = d$tau[index_max],
        setting_id_at_max = as.integer(d$setting_id[index_max]),
        theta_assumed_at_max = d$theta_assumed[index_max],
        R_global_unrounded = as.integer(global_required),
        R_global_rounded_100 = as.integer(
          ceiling(global_required / 100) * 100
        )
      ))
    }
  ))

  valid_R <- setting_summary$R_required_sensitivity[
    is.finite(setting_summary$R_required_sensitivity)
  ]

  # Step 5: select the largest setting-specific requirement and round upward.
  # Upward rounding cannot reduce the planned precision.
  R_global <- if (length(valid_R) > 0L) max(valid_R) else NA_integer_
  R_rounded_100 <- if (is.finite(R_global)) {
    as.integer(ceiling(R_global / 100) * 100)
  } else {
    NA_integer_
  }

  # Step 6: store the assumptions and final choice in a compact table suitable
  # for direct inspection or inclusion in a report appendix.
  global_summary <- tibble(
    design_name = design$design_name,
    n_settings = nrow(settings),
    n_tasks = nrow(run_results),
    n_success = sum(run_results$ok),
    n_failed = sum(!run_results$ok),
    mcse_target = design$mcse_target,
    R_min = design$R_min,
    sd_multiplier_min = min(sd_multiplier_grid),
    sd_multiplier_max = max(sd_multiplier_grid),
    sd_multiplier_grid = paste(sd_multiplier_grid, collapse = ","),
    theta_assumed_min = min(theta_assumed_grid),
    theta_assumed_max = max(theta_assumed_grid),
    max_raw_tau_rf_lasso = safe_max(setting_summary$raw_tau_rf_vs_lasso),
    max_raw_R_required = safe_max(setting_summary$raw_R_required),
    max_sensitivity_tau = safe_max(setting_summary$surrogate_tau_max),
    R_global_sensitivity = R_global,
    R_global_rounded_100 = R_rounded_100,
    n_test_policy = design$n_test_policy,
    note = paste(
      "R_global_sensitivity uses the Random Forest correlation pattern and",
      "fixed comparator standard-deviation scenarios relative to the LASSO",
      "standard deviation. It is not an empirical Ridge-to-LASSO",
      "standard-deviation estimate."
    ),
    created_at = Sys.time()
  )

  return(list(
    setting_summary = setting_summary,
    sensitivity_summary = sensitivity_summary,
    global_summary = global_summary
  ))
}

###############################################################
# E. COMPLETE STANDARD-DEVIATION AND PRECISION PILOT
###############################################################
# PURPOSE
# Run the complete blinded precision pilot from design construction through
# saved tables and the report figure.
#
# WHY THIS FUNCTION EXISTS
# This is the orchestration function used by run/02_run_precision_pilot.R.
# It connects the general simulation functions in pilot_core.R with the
# precision-specific summary above while keeping execution choices explicit.
#
# INPUTS
# R_sd: Number of pilot repetitions per simulation setting.
# output_dir: Destination for final output files and, by default, the registry.
# backend: Use "batchtools" for resumable parallel execution. Use "serial" for
#   a simple sequential run.
# n_cores: Number of local workers used by the batchtools backend.
# seed_base: Starting value used to create deterministic task seeds.
# design: Shared simulation constants returned by make_design().
# registry_dir: Persistent batchtools registry used for resuming a run.
# resume, overwrite_registry, resubmit_errors, wait: Execution controls passed
#   to run_tasks(), their behavior is documented in pilot_core.R.
# save_outputs: If TRUE, save a summary R data object, three comma-separated
#   summary tables, and the figure in Portable Network Graphics and Portable
#   Document Format. Repetition-level prediction errors are excluded.
# sd_multiplier_grid, theta_assumed_grid: Precision sensitivity assumptions.
# source_files: R files that batchtools workers must source before they execute
#   a task. The worker function is defined in pilot_core.R, so the launcher
#   correctly supplies only that core file.
#
# RETURN VALUE
# A list containing run_stamp, the design, settings, execution status,
# summaries, sensitivity assumptions, and, when save_outputs is TRUE, the paths
# of all written files. The returned object deliberately excludes
# repetition-level prediction errors. run_stamp is NA when save_outputs is FALSE
# because no persistent run file family is created.
run_precision_pilot <- function(
    R_sd = 500L,
    output_dir = file.path("outputs", "precision_pilot"),
    backend = c("batchtools", "serial"),
    n_cores = 8L,
    seed_base = 1000L,
    design = make_design(),
    registry_dir = file.path(output_dir, "registry"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = FALSE,
    wait = TRUE,
    save_outputs = TRUE,
    sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3),
    theta_assumed_grid = c(-0.02, 0, 0.02),
    source_files = script_file) {

  # Step 1: resolve the execution mode and check the plotting dependency.
  backend <- match.arg(backend)

  # Check ggplot2 before starting a pilot with saved figures.
  if (isTRUE(save_outputs) && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("save_outputs = TRUE requires package 'ggplot2'.")
  }

  # Step 2: expand the fixed factorial design and create one deterministic task
  # for every setting-by-repetition combination.
  settings <- make_all_settings(design)
  tasks <- make_tasks(settings, R = R_sd, seed_base = seed_base)

  # Step 3: execute or resume the task plan. Each task fits LASSO and Random
  # Forest to the same simulated training/test data.
  task_run <- run_tasks(
    tasks = tasks,
    task_function = run_one_sd_task,
    design = design,
    backend = backend,
    registry_dir = registry_dir,
    n_cores = n_cores,
    cluster_function = design$batchtools_cluster_function,
    resume = resume,
    overwrite_registry = overwrite_registry,
    resubmit_errors = resubmit_errors,
    wait = wait,
    source_files = source_files
  )

  # Step 4: convert paired task results into setting-specific and global
  # precision-planning summaries.
  summaries <- summarise_precision_pilot(
    run_results = task_run$run_results,
    settings = settings,
    design = design,
    sd_multiplier_grid = sd_multiplier_grid,
    theta_assumed_grid = theta_assumed_grid
  )

  # Step 5: assemble the public summary object. The paired prediction errors
  # have already been converted into standard deviations, coefficients of
  # variation, dependence measures, influence-value standard deviations, and
  # repetition requirements. They are therefore omitted from both the returned
  # object and the final R data file. During a batchtools run, the local
  # registry still holds task results so an interrupted calculation can resume.
  # That registry is execution state and is not part of the public deposit.
  out <- list(
    run_stamp = NA_character_,
    design = design,
    settings = settings,
    backend = backend,
    seed_base = as.integer(seed_base),
    batchtools_status = task_run$status,
    setting_summary = summaries$setting_summary,
    sensitivity_summary = summaries$sensitivity_summary,
    global_summary = summaries$global_summary,
    precision_sensitivity = list(
      sd_multiplier_grid = sd_multiplier_grid,
      theta_assumed_grid = theta_assumed_grid
    )
  )

  # Step 6: save timestamped numerical results and a timestamped figure. The
  # figure from this execution includes its timestamp. The figure-only
  # launcher accepts the exact setting-summary path independently.
  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- make_timestamp()
    f_rds <- file.path(output_dir, paste0("sd_precision_pilot_", stamp, ".rds"))
    f_set <- file.path(output_dir, paste0("sd_precision_setting_summary_", stamp, ".csv"))
    f_sensitivity <- file.path(
      output_dir,
      paste0("sd_precision_sensitivity_summary_", stamp, ".csv")
    )
    f_global <- file.path(output_dir, paste0("sd_precision_global_summary_", stamp, ".csv"))
    figure_stem <- paste0(
      "Supp_Fig_RF_LASSO_SD_ratio_boxplot_pretty_",
      stamp
    )

    figure <- save_rf_lasso_sd_ratio_figure(
      out$setting_summary,
      output_dir = output_dir,
      file_stem = figure_stem
    )
    out$run_stamp <- stamp
    out$rf_lasso_sd_figure_summary <- figure$summary
    out$saved_files <- c(
      rds = f_rds,
      setting_summary = f_set,
      sensitivity_summary = f_sensitivity,
      global_summary = f_global,
      figure_png = figure$png,
      figure_pdf = figure$pdf
    )

    # Record the file paths inside the saved object as well as in the returned
    # object. Write the readable tables before the R data file so the summary
    # object is created only after its referenced numerical extracts succeed.
    write.csv(out$setting_summary, f_set, row.names = FALSE)
    write.csv(out$sensitivity_summary, f_sensitivity, row.names = FALSE)
    write.csv(out$global_summary, f_global, row.names = FALSE)
    saveRDS(out, f_rds)
  }

  return(out)
}

###############################################################
# EXAMPLE: COMPLETE PILOT
###############################################################
# The recommended user-facing entry point is run/02_run_precision_pilot.R. The
# direct function call below is useful only when custom execution arguments
# are required. Run package installation and tests before this call.
#
# sd_pilot <- run_precision_pilot(
#   R_sd = 500L,
#   output_dir = file.path("outputs", "precision_pilot"),
#   backend = "batchtools",
#   n_cores = 8L,
#   registry_dir = file.path(
#     "outputs", "precision_pilot", "registry"
#   ),
#   save_outputs = TRUE,
#   source_files = script_file
# )
###############################################################
