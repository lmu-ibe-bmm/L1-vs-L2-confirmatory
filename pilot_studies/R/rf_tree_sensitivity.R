# Check whether num.trees changes the RF/LASSO variability used for precision planning.

# Locate a shared module from the project or the active source-file directory.
# Inputs: filename. Returns an existing absolute path or stops.
find_rf_sensitivity_sibling <- function(filename) {

  candidates <- file.path("R", filename)

  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("--file=", command_line, fixed = TRUE)
  if (length(file_argument) > 0L) {
    current_file <- sub(
      "--file=", "", command_line[file_argument[1]], fixed = TRUE
    )
    current_file <- normalizePath(
      current_file, winslash = "/", mustWork = FALSE
    )
    candidates <- c(candidates, file.path(dirname(current_file), filename))
  }

  source_frames <- sys.frames()
  source_files <- vapply(
    source_frames,

    function(frame) {
      return(if (!is.null(frame$ofile)) as.character(frame$ofile) else NA_character_)
    },
    character(1)
  )
  source_files <- source_files[!is.na(source_files)]
  if (length(source_files) > 0L) {
    source_files <- vapply(
      source_files,
      normalizePath,
      character(1),
      winslash = "/",
      mustWork = FALSE
    )

    module_files <- source_files[
      basename(source_files) == "rf_tree_sensitivity.R"
    ]
    candidates <- c(
      if (length(module_files) > 0L) {
        file.path(dirname(module_files), filename)
      } else {
        character(0)
      },
      candidates,
      file.path(dirname(source_files), filename)
    )
  }

  candidates <- unique(c(candidates, filename))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("File not found: ", filename)
  }

  return(normalizePath(existing[1], winslash = "/", mustWork = TRUE))
}

rf_tree_sensitivity_script_file <- find_rf_sensitivity_sibling(
  "rf_tree_sensitivity.R"
)
rf_tree_sensitivity_precision_file <- find_rf_sensitivity_sibling(
  "precision_pilot.R"
)
rf_tree_sensitivity_core_file <- find_rf_sensitivity_sibling(
  "pilot_core.R"
)

if (!exists("summarise_precision_pilot", mode = "function")) {
  source(rf_tree_sensitivity_precision_file)
}

# Fit LASSO once and compare RF num.trees values within the same simulated data.
# Inputs: one task_row and design with rf_tree_counts.
# Returns one row with paired test MSEs, fit times and any failure details.
# Data, LASSO folds and RF initialization seed are held fixed across tree counts.
run_one_rf_tree_sensitivity_task <- function(
    task_row,
    design = make_design()) {

  task_start <- proc.time()[["elapsed"]]
  task_row <- tibble::as_tibble(task_row)

  if (is.null(design$rf_tree_counts)) {
    stop("design$rf_tree_counts is required for this sensitivity analysis.")
  }

  tree_counts <- sort(unique(as.integer(design$rf_tree_counts)))
  if (
    length(tree_counts) < 2L ||
      any(!is.finite(tree_counts)) ||
      any(tree_counts < 1L)
  ) {
    stop("rf_tree_counts must contain at least two positive integers.")
  }

  result <- base_result_template(task_row) %>%
    mutate(
      mse_lasso = NA_real_,
      lambda_lasso_min = NA_real_,
      mtry = NA_integer_,
      time_lasso_fit_sec = NA_real_,
      time_lasso_pred_sec = NA_real_
    )

  for (tree_count in tree_counts) {
    result[[paste0("mse_rf_", tree_count)]] <- NA_real_
    result[[paste0("time_rf_fit_", tree_count, "_sec")]] <- NA_real_
    result[[paste0("time_rf_pred_", tree_count, "_sec")]] <- NA_real_
  }

  active_step <- "initialization"

  result <- tryCatch({

    active_step <- "beta"
    beta_start <- proc.time()[["elapsed"]]
    beta_object <- make_beta(
      task_row,
      seed = task_row$seed + 1L,
      design = design
    )
    result$time_beta_sec <- proc.time()[["elapsed"]] - beta_start
    result <- fill_beta_info(result, beta_object)

    active_step <- "data"
    simulated_data <- simulate_data(
      task_row,
      beta_object = beta_object,
      seed = task_row$seed + 2L,
      design = design
    )
    result$time_data_sec <- simulated_data$time_data_sec

    active_step <- "foldid"
    foldid <- make_cv_foldid(
      n_train = nrow(simulated_data$X_train),
      nfolds = design$glmnet_nfolds,
      seed = task_row$seed + 3L
    )

    active_step <- "lasso"
    lasso_fit <- fit_glmnet(
      simulated_data,
      alpha = 1,
      foldid = foldid,
      beta_true = beta_object$beta,
      active_variables = beta_object$active_variables,
      design = design
    )
    result$mse_lasso <- lasso_fit$mse
    result$lambda_lasso_min <- lasso_fit$lambda_min
    result$time_lasso_fit_sec <- lasso_fit$fit_time_sec
    result$time_lasso_pred_sec <- lasso_fit$pred_time_sec

    predictor_count <- as.integer(task_row$p)
    mtry_value <- max(1L, as.integer(floor(sqrt(predictor_count))))
    result$mtry <- mtry_value

    for (tree_count in tree_counts) {
      active_step <- paste0("rf_fit_", tree_count, "_trees")
      fit_start <- proc.time()[["elapsed"]]
      forest_fit <- ranger::ranger(
        x = simulated_data$X_train,
        y = simulated_data$y_train,
        num.trees = as.integer(tree_count),
        mtry = mtry_value,
        min.node.size = as.integer(design$rf_min_node_size),
        seed = as.integer(task_row$seed + 4L),
        num.threads = 1L,
        write.forest = TRUE
      )
      result[[paste0("time_rf_fit_", tree_count, "_sec")]] <-
        proc.time()[["elapsed"]] - fit_start

      active_step <- paste0("rf_prediction_", tree_count, "_trees")
      prediction_start <- proc.time()[["elapsed"]]
      forest_prediction <- as.numeric(stats::predict(
        forest_fit,
        data = simulated_data$X_test,
        num.threads = 1L
      )$predictions)
      result[[paste0("time_rf_pred_", tree_count, "_sec")]] <-
        proc.time()[["elapsed"]] - prediction_start

      performance <- compute_prediction_metrics(
        simulated_data$y_test,
        forest_prediction
      )
      result[[paste0("mse_rf_", tree_count)]] <- performance$mse

      rm(forest_fit, forest_prediction)
    }

    pmse_columns <- paste0("mse_rf_", tree_counts)
    if (
      !is.finite(result$mse_lasso) ||
        any(!is.finite(unlist(result[pmse_columns], use.names = FALSE)))
    ) {
      stop("At least one prediction mean squared error is not finite.")
    }
    result$ok <- TRUE
    result

  }, error = function(error_condition) {

    result$failed_step <- active_step
    result$error_message <- conditionMessage(error_condition)
    return(result)
  })

  result$time_total_sec <- proc.time()[["elapsed"]] - task_start
  return(result)
}

# Apply the precision-planning formula separately to each RF tree count.
# Inputs: run_results, settings, design and the SD/effect scenario grids.
# Returns setting/global summaries and within-setting changes from 200 trees.
# Planning uses hypothetical comparator SDs and RF-LASSO correlation, not Ridge fits.
summarise_rf_tree_sensitivity <- function(
    run_results,
    settings,
    design = make_design(),
    sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3),
    theta_assumed_grid = c(-0.02, 0, 0.02)) {

  tree_counts <- sort(unique(as.integer(design$rf_tree_counts)))
  baseline_trees <- 200L

  required_columns <- c(
    "task_id", "setting_id", "run_id", "ok", "mse_lasso",
    paste0("mse_rf_", tree_counts)
  )
  missing_columns <- setdiff(required_columns, names(run_results))
  if (length(missing_columns) > 0L) {
    stop(
      "Tree-count sensitivity results are missing required columns: ",
      paste(missing_columns, collapse = ", "), "."
    )
  }
  if (any(!run_results$ok)) {
    stop(
      "The tree-count sensitivity run contains failed tasks and cannot be ",
      "summarized."
    )
  }
  if (anyDuplicated(run_results$task_id) > 0L) {
    stop("The tree-count sensitivity results contain duplicate task IDs.")
  }
  if (!baseline_trees %in% tree_counts) {
    stop("The sensitivity analysis must include the 200-tree baseline.")
  }

  setting_tables <- vector("list", length(tree_counts))
  sensitivity_tables <- vector("list", length(tree_counts))
  global_tables <- vector("list", length(tree_counts))

  for (tree_index in seq_along(tree_counts)) {
    tree_count <- tree_counts[tree_index]
    forest_column <- paste0("mse_rf_", tree_count)

    precision_input <- run_results[, c(
      "task_id", "setting_id", "run_id", "ok", "mse_lasso", forest_column
    ), drop = FALSE]
    names(precision_input)[names(precision_input) == forest_column] <- "mse_rf"

    precision_summary <- summarise_precision_pilot(
      run_results = precision_input,
      settings = settings,
      design = design,
      sd_multiplier_grid = sd_multiplier_grid,
      theta_assumed_grid = theta_assumed_grid
    )

    setting_table <- precision_summary$setting_summary
    setting_table$num_trees <- as.integer(tree_count)
    setting_table$sd_ratio_rf_lasso <-
      setting_table$sd_mse_rf_raw / setting_table$sd_mse_lasso
    setting_tables[[tree_index]] <- setting_table %>%
      relocate(num_trees, .after = setting_id)

    sensitivity_table <- precision_summary$sensitivity_summary
    sensitivity_table$num_trees <- as.integer(tree_count)
    sensitivity_tables[[tree_index]] <- sensitivity_table %>%
      relocate(num_trees)

    global_table <- precision_summary$global_summary
    global_table$num_trees <- as.integer(tree_count)
    global_tables[[tree_index]] <- global_table %>%
      relocate(num_trees)
  }

  setting_summary <- bind_rows(setting_tables)
  sensitivity_summary <- bind_rows(sensitivity_tables)
  global_summary <- bind_rows(global_tables)

  tree_count_summary <- setting_summary %>%
    group_by(num_trees) %>%
    summarise(
      n_settings = n(),
      median_sd_ratio_rf_lasso = stats::median(sd_ratio_rf_lasso),
      q25_sd_ratio_rf_lasso = as.numeric(stats::quantile(
        sd_ratio_rf_lasso, probs = 0.25, names = FALSE
      )),
      q75_sd_ratio_rf_lasso = as.numeric(stats::quantile(
        sd_ratio_rf_lasso, probs = 0.75, names = FALSE
      )),
      max_sd_ratio_rf_lasso = max(sd_ratio_rf_lasso),
      n_ratio_above_1 = sum(sd_ratio_rf_lasso > 1),
      n_ratio_above_1_3 = sum(sd_ratio_rf_lasso > 1.3),
      n_ratio_above_2 = sum(sd_ratio_rf_lasso > 2),
      min_correlation_lasso_rf = min(corr_lasso_rf_raw),
      median_correlation_lasso_rf = stats::median(corr_lasso_rf_raw),
      max_correlation_lasso_rf = max(corr_lasso_rf_raw),
      .groups = "drop"
    ) %>%
    left_join(
      global_summary %>%
        select(
          num_trees,
          max_sensitivity_tau,
          R_global_sensitivity,
          R_global_rounded_100
        ),
      by = "num_trees"
    )

  baseline <- setting_summary %>%
    filter(num_trees == baseline_trees) %>%
    select(
      setting_id,
      baseline_sd_mse_rf = sd_mse_rf_raw,
      baseline_sd_ratio_rf_lasso = sd_ratio_rf_lasso,
      baseline_corr_lasso_rf = corr_lasso_rf_raw,
      baseline_raw_R_required = raw_R_required,
      baseline_R_required_sensitivity = R_required_sensitivity
    )

  comparison_with_200 <- setting_summary %>%
    select(
      setting_id, num_trees,
      n, p_over_n, p, rho, a, k, snr,
      sd_mse_lasso, sd_mse_rf_raw, sd_ratio_rf_lasso,
      corr_lasso_rf_raw, raw_R_required, R_required_sensitivity
    ) %>%
    left_join(baseline, by = "setting_id") %>%
    mutate(
      change_sd_mse_rf_from_200 =
        sd_mse_rf_raw - baseline_sd_mse_rf,
      change_sd_ratio_from_200 =
        sd_ratio_rf_lasso - baseline_sd_ratio_rf_lasso,
      change_correlation_from_200 =
        corr_lasso_rf_raw - baseline_corr_lasso_rf,
      change_raw_R_from_200 =
        raw_R_required - baseline_raw_R_required,
      change_sensitivity_R_from_200 =
        R_required_sensitivity - baseline_R_required_sensitivity
    )

  return(list(
    setting_summary = setting_summary,
    tree_count_summary = tree_count_summary,
    sensitivity_summary = sensitivity_summary,
    global_summary = global_summary,
    comparison_with_200 = comparison_with_200
  ))
}

# Plot the setting-level RF/LASSO SD ratios for each num.trees condition.
# Inputs: setting_summary, output_dir, file_stem and image options.
# Returns PNG/PDF paths and the plot. Lines at 1 and 1.3 mark planning scenarios.
save_rf_tree_sensitivity_figure <- function(
    setting_summary,
    output_dir,
    file_stem,
    width_in = 8.0,
    height_in = 5.5,
    dpi = 300) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for the sensitivity figure.")
  }

  required_columns <- c(
    "setting_id", "num_trees", "sd_ratio_rf_lasso"
  )
  missing_columns <- setdiff(required_columns, names(setting_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "The figure input is missing required columns: ",
      paste(missing_columns, collapse = ", "), "."
    )
  }

  plot_data <- setting_summary[, required_columns, drop = FALSE]
  if (any(!is.finite(plot_data$sd_ratio_rf_lasso))) {
    stop("All standard-deviation ratios must be finite before plotting.")
  }

  tree_levels <- sort(unique(as.integer(plot_data$num_trees)))
  plot_data$tree_count <- factor(
    plot_data$num_trees,
    levels = tree_levels,
    labels = format(
      tree_levels, big.mark = ",", scientific = FALSE, trim = TRUE
    )
  )

  sensitivity_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = tree_count, y = sd_ratio_rf_lasso)
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      colour = "#1F77B4",
      linewidth = 0.75,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 1.3,
      colour = "#1F77B4",
      linewidth = 0.75,
      linetype = "dotted"
    ) +
    ggplot2::geom_boxplot(
      width = 0.55,
      fill = "white",
      colour = "black",
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.10, height = 0),
      shape = 21,
      size = 1.5,
      stroke = 0.35,
      fill = "white",
      colour = "black",
      alpha = 0.42
    ) +
    ggplot2::labs(
      x = "Number of trees",
      y = "Ratio of setting-specific SD(PMSE): Random Forest / LASSO"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.7
      ),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_file <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_stem, ".pdf"))

  ggplot2::ggsave(
    png_file,
    plot = sensitivity_plot,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  ggplot2::ggsave(
    pdf_file,
    plot = sensitivity_plot,
    width = width_in,
    height = height_in,
    units = "in",
    device = pdf_device,
    bg = "white"
  )

  return(list(
    png = png_file,
    pdf = pdf_file,
    plot = sensitivity_plot
  ))
}

# Run or resume the paired RF sensitivity check with corrected component seeds.
# Inputs: repetitions, tree_counts, optional setting_ids and execution/save options.
# Returns summary tables and provenance. The full launcher uses all 192 settings.
# save_outputs writes the RDS object, CSV tables and SD-ratio figure.
run_rf_tree_sensitivity <- function(
    R_sensitivity = 100L,
    tree_counts = c(200L, 500L, 1000L),
    setting_ids = NULL,
    output_dir = file.path("outputs", "rf_tree_sensitivity"),
    backend = c("batchtools", "serial"),
    n_cores = 8L,
    seed_base = 400000000L,
    design = make_design(),
    registry_dir = file.path(output_dir, "registry"),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = FALSE,
    wait = TRUE,
    save_outputs = TRUE,
    sd_multiplier_grid = c(1.0, 1.1, 1.2, 1.3),
    theta_assumed_grid = c(-0.02, 0, 0.02),
    source_files = rf_tree_sensitivity_core_file,
    seed_policy = "separated") {

  backend <- match.arg(backend)
  seed_policy <- match.arg(seed_policy, "separated")
  R_sensitivity <- as.integer(R_sensitivity)
  tree_counts <- sort(unique(as.integer(tree_counts)))

  if (!is.finite(R_sensitivity) || R_sensitivity < 2L) {
    stop("R_sensitivity must be an integer of at least two.")
  }
  if (
    length(tree_counts) < 2L ||
      any(!is.finite(tree_counts)) ||
      any(tree_counts < 1L)
  ) {
    stop("tree_counts must contain at least two positive integers.")
  }
  if (!200L %in% tree_counts) {
    stop("tree_counts must include the original 200-tree condition.")
  }
  if (isTRUE(save_outputs) && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("save_outputs = TRUE requires package 'ggplot2'.")
  }

  design$rf_tree_counts <- tree_counts
  design$rf_tree_sensitivity_baseline <- 200L
  design$rf_tree_sensitivity_seed_policy <- seed_policy

  design$rf_tree_sensitivity_code_md5 <- unname(tools::md5sum(c(
    rf_tree_sensitivity_core_file,
    rf_tree_sensitivity_precision_file,
    rf_tree_sensitivity_script_file
  )))
  if (anyNA(design$rf_tree_sensitivity_code_md5)) {
    stop("All sensitivity source files must be readable before execution.")
  }

  all_settings <- make_all_settings(design)
  if (is.null(setting_ids)) {
    settings <- all_settings
    setting_selection <- "all_settings"
  } else {
    setting_ids <- sort(unique(as.integer(setting_ids)))
    if (
      length(setting_ids) < 1L ||
        any(!setting_ids %in% all_settings$setting_id)
    ) {
      stop("Every setting_ids value must identify a design setting.")
    }
    settings <- all_settings[
      all_settings$setting_id %in% setting_ids,
      ,
      drop = FALSE
    ]
    setting_selection <- "selected_settings"
  }

  tasks <- make_tasks(settings = settings, R = R_sensitivity, seed_base = seed_base)

  task_run <- run_tasks(
    tasks = tasks,
    task_function = run_one_rf_tree_sensitivity_task,
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

  summaries <- summarise_rf_tree_sensitivity(
    run_results = task_run$run_results,
    settings = settings,
    design = design,
    sd_multiplier_grid = sd_multiplier_grid,
    theta_assumed_grid = theta_assumed_grid
  )

  public_result <- list(
    run_stamp = NA_character_,
    design = design,
    settings = settings,
    setting_selection = setting_selection,
    R_sensitivity = R_sensitivity,
    tree_counts = tree_counts,
    backend = backend,
    seed_base = as.integer(seed_base),
    batchtools_status = task_run$status,
    setting_summary = summaries$setting_summary,
    tree_count_summary = summaries$tree_count_summary,
    sensitivity_summary = summaries$sensitivity_summary,
    global_summary = summaries$global_summary,
    comparison_with_200 = summaries$comparison_with_200,
    precision_sensitivity = list(
      sd_multiplier_grid = sd_multiplier_grid,
      theta_assumed_grid = theta_assumed_grid
    )
  )

  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    run_stamp <- make_timestamp()

    rds_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_", run_stamp, ".rds")
    )
    setting_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_setting_summary_", run_stamp, ".csv")
    )
    tree_count_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_tree_count_summary_", run_stamp, ".csv")
    )
    sensitivity_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_precision_summary_", run_stamp, ".csv")
    )
    global_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_global_summary_", run_stamp, ".csv")
    )
    comparison_file <- file.path(
      output_dir,
      paste0("rf_tree_sensitivity_comparison_with_200_", run_stamp, ".csv")
    )

    figure_result <- save_rf_tree_sensitivity_figure(
      setting_summary = public_result$setting_summary,
      output_dir = output_dir,
      file_stem = paste0("rf_tree_sensitivity_figure_", run_stamp)
    )

    public_result$run_stamp <- run_stamp
    public_result$saved_files <- c(
      rds = rds_file,
      setting_summary = setting_file,
      tree_count_summary = tree_count_file,
      precision_summary = sensitivity_file,
      global_summary = global_file,
      comparison_with_200 = comparison_file,
      figure_png = figure_result$png,
      figure_pdf = figure_result$pdf
    )

    utils::write.csv(
      public_result$setting_summary,
      setting_file,
      row.names = FALSE
    )
    utils::write.csv(
      public_result$tree_count_summary,
      tree_count_file,
      row.names = FALSE
    )
    utils::write.csv(
      public_result$sensitivity_summary,
      sensitivity_file,
      row.names = FALSE
    )
    utils::write.csv(
      public_result$global_summary,
      global_file,
      row.names = FALSE
    )
    utils::write.csv(
      public_result$comparison_with_200,
      comparison_file,
      row.names = FALSE
    )
    saveRDS(public_result, rds_file)
  }

  return(public_result)
}
