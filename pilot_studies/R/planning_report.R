# Read complete pilots and update the planning and runtime tables.

source("R/workflow_helpers.R")
source("R/precision_pilot.R")

# Read a completed precision/RF result after checking its source and output hashes.
# Inputs: output_dir. Loads and validates the saved summary files.
read_completed_run <- function(output_dir) {
  file <- file.path(output_dir, "run_complete.rds")
  if (!file.exists(file)) stop("This run is not complete yet: ", output_dir)
  record <- readRDS(file)
  result <- load_completed_run(output_dir, record$specification,
                                 names(record$source_md5))
  if (is.null(result)) stop("Saved result files are missing. Resume the run first.")
  return(result)
}

# Select the latest compatible runtime benchmark, preferring outputs over retained results.
# Inputs: none. Returns the benchmark with its source-file path.
read_runtime_pilot <- function() {
  paths <- c("outputs/runtime_pilot", "results/runtime_pilot")
  for (path in paths) {
    files <- list.files(path, pattern = "^runtime_pilot_.*[.]rds$", full.names = TRUE)
    if (!length(files)) next
    file <- tail(sort(files), 1L)
    result <- readRDS(file)
    if (nrow(result$settings) != 192L ||
        !isTRUE(all.equal(result$design$en_alpha, 1 / 3)) ||
        !isTRUE(all.equal(result$design$en_near_lasso_alpha, 0.95)) ||
        result$design$n_test != 5000L || result$runtime_global$n_failed != 0L ||
        anyNA(result$runtime_by_setting$mean_time_total_sec) ||
        any(result$runtime_by_setting$n_success != result$runtime_global$R_runtime)) {
      stop("The runtime benchmark does not match the corrected setup.")
    }
    result$source_file <- file
    return(result)
  }
  stop("No corrected runtime benchmark is available.")
}

# Compare the pilot requirement and rounding options with the configured study R.
# Inputs: the complete precision result and optional main_R.
# Returns candidate R values and tau / sqrt(R). Rounding is not an uncertainty bound.
make_planning_candidates <- function(precision, main_R = NA_integer_) {
  check_precision_run(precision)
  point <- max(precision$setting_summary$R_required_sensitivity)
  tau <- precision$global_summary$max_sensitivity_tau
  if (length(main_R) != 1L || (!is.na(main_R) &&
      (!is.finite(main_R) || main_R < 1 || main_R != floor(main_R)))) {
    stop("main_R must be NA or one positive integer.")
  }
  result <- data.frame(
    scenario = c("pilot_requirement", "round_up_100", "round_up_500",
                 "round_up_1000", "comparison_5200", "comparison_5800"),
    repetitions_per_setting = as.integer(c(point, ceiling(point / 100) * 100,
       ceiling(point / 500) * 500, ceiling(point / 1000) * 1000, 5200, 5800)),
    stringsAsFactors = FALSE
  )
  if (!is.na(main_R)) result <- rbind(result, data.frame(
    scenario = "selected_main_R", repetitions_per_setting = as.integer(main_R)))
  result$estimated_planning_mcse <- tau / sqrt(result$repetitions_per_setting)
  result$mcse_target <- precision$design$mcse_target
  result$meets_pilot_point_target <- result$estimated_planning_mcse <= result$mcse_target
  result$is_selected_main_R <- result$scenario == "selected_main_R"
  return(result)
}

# Combine completed pilots with the runtime benchmark without fitting models.
# Inputs: output_dir, with main_R and worker count read from project_config.
# Writes the selected plan, comparison tables and status. Returns candidates invisibly.
summarize_planning <- function(output_dir = "outputs/planning_summary") {
  precision <- read_completed_run("outputs/precision_pilot_R500")
  rf <- read_completed_run("outputs/rf_tree_sensitivity_R100")
  check_precision_run(precision)
  check_rf_run(rf)
  if (precision$seed_base != project_config$precision_seed_base ||
      rf$seed_base != project_config$rf_seed_base) {
    stop("The completed runs do not use the recorded seed families.")
  }
  runtime <- read_runtime_pilot()
  columns <- c("setting_id", "n", "p_over_n", "p", "rho", "a", "k", "snr", "b", "n_test")
  if (!isTRUE(all.equal(as.data.frame(runtime$settings[, columns]),
                       as.data.frame(precision$settings[, columns]),
                       check.attributes = FALSE))) {
    stop("Precision and runtime grids do not match.")
  }
  candidates <- make_planning_candidates(precision, project_config$main_R)
  times <- dplyr::bind_rows(lapply(candidates$repetitions_per_setting, function(R) {
    project_planned_study_runtime(runtime$runtime_by_setting, runtime$runtime_global,
                                 repetitions_per_setting = R, n_cores = project_config$n_cores)
  }))
  candidates <- cbind(candidates, times[, setdiff(names(times), "repetitions_per_setting")])
  candidates$runtime_source <- runtime$source_file
  candidates$precision_source <- precision$saved_files[["rds"]]
  candidates$runtime_days_mean <- "24 hours of continuous operation at the stated worker count"
  selected <- candidates[candidates$is_selected_main_R, , drop = FALSE]
  if (nrow(selected)) {
    selected$selection_rule <- "Maximum over all 192 settings in the R500 pilot, rounded upward to 100"
    selected$pilot_repetitions_per_setting <- 500L
    selected$pilot_maximum_integer_requirement <- max(precision$setting_summary$R_required_sensitivity)
    selected$matches_round_up_100 <- selected$repetitions_per_setting ==
      ceiling(selected$pilot_maximum_integer_requirement / 100) * 100
  }
  top <- precision$setting_summary[
    order(-precision$setting_summary$R_required_sensitivity,
          precision$setting_summary$setting_id), ]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(precision$global_summary,
                  file.path(output_dir, "precision_global_summary.csv"), row.names = FALSE)
  utils::write.csv(head(top, 20L),
                  file.path(output_dir, "precision_top_settings.csv"), row.names = FALSE)
  utils::write.csv(candidates,
                  file.path(output_dir, "precision_runtime_candidates.csv"), row.names = FALSE)
  utils::write.csv(selected,
                  file.path(output_dir, "main_study_plan.csv"), row.names = FALSE)
  utils::write.csv(rf$tree_count_summary,
                  file.path(output_dir, "rf_tree_count_summary.csv"), row.names = FALSE)
  point <- max(top$R_required_sensitivity)
  text <- c(
    "192 settings, 500 completed precision repetitions per setting.",
    paste("Maximum integer R requirement:", point),
    paste("Maximum rounded upward to 100:", ceiling(point / 100) * 100),
    paste("Setting ID at the maximum:", top$setting_id[1L]),
    if (is.na(project_config$main_R)) "The main-study R is still open." else
      paste("Selected main-study R:", project_config$main_R),
    "Selection rule: complete R500 pilot maximum, rounded upward to 100.",
    "The separate R100 RF check is a sensitivity analysis, not the selection pilot.",
    "Rounding scenarios are not uncertainty bounds or automatic decisions.",
    "These are RF/LASSO-based planning values, not measured Ridge-LASSO MCSEs.",
    "Runtime days assume continuous operation. Shutdowns add calendar time.",
    paste("Runtime source:", runtime$source_file)
  )
  writeLines(text, file.path(output_dir, "planning_status.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  cat(paste(text, collapse = "\n"), "\n")
  print(candidates[, c("scenario", "repetitions_per_setting", "estimated_planning_mcse",
                       "projected_ideal_parallel_days", "projected_adjusted_days")])
  return(invisible(candidates))
}
