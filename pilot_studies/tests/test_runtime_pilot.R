###############################################################
# TESTS FOR THE RUNTIME PILOT
###############################################################
# PURPOSE
# These tests verify one complete timed task with all four glmnet models and
# the algebra used to project pilot timings to a larger simulation. They focus
# on measured components and projection logic, not on rechecking the fixed
# number of design settings or tasks.
#
# WHY THESE TESTS MATTER
# Runtime projections can be wrong even when every model fit succeeds. The
# smoke test confirms that one task records each required timing component.
# The synthetic examples then make the projection arithmetic transparent and
# verify that scheduler overhead is used only when the timestamp information
# needed to estimate it is complete.
#
# WHEN TO RUN
# Run through tests/run_all_tests.R after package installation and before
# run/03_run_runtime_pilot.R. To run only this file from the top-level package
# directory:
#   Rscript tests/test_runtime_pilot.R
# The model-fitting smoke test uses a reduced test set and does not launch
# batchtools or the full runtime pilot.
###############################################################

source(file.path("R", "runtime_pilot.R"))

# Reduce only the independent test-set size for the smoke test. The demanding
# correlation and support structure remain unchanged.
smoke_design <- design
smoke_design$n_test <- 200L

smoke_setting <- make_setting(
  n = 200,
  p_over_n = 0.2,
  rho = 0.95,
  a = 0.75,
  k = 10,
  snr = 2,
  design = smoke_design
)
smoke_setting$setting_id <- 1L
smoke_task <- make_tasks(smoke_setting, R = 1L, seed_base = 9000L)

packages_available <- all(vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
))

if (packages_available) {
  task_result <- run_one_runtime_task(smoke_task, smoke_design)
  stopifnot(
    isTRUE(task_result$ok),
    all(is.finite(c(
      task_result$time_total_sec,
      task_result$time_data_sec,
      task_result$time_X_train_sec,
      task_result$time_X_test_sec,
      task_result$time_ridge_fit_sec,
      task_result$time_lasso_fit_sec,
      task_result$time_en_fit_sec,
      task_result$time_en095_fit_sec,
      task_result$time_prediction_metrics_sec,
      task_result$time_coefficient_metrics_sec,
      task_result$lambda_ridge_min,
      task_result$lambda_lasso_min,
      task_result$lambda_en_min,
      task_result$lambda_en095_min
    ))),
    task_result$time_total_sec >= 0,
    task_result$time_internal_unaccounted_sec >= 0
  )

  # Summarize internal task times with an empty scheduler table.
  runtime_summary <- summarise_runtime_detail(
    run_results = task_result,
    settings = smoke_setting,
    job_detail = tibble(
      submitted = as.POSIXct(character()),
      started = as.POSIXct(character()),
      done = as.POSIXct(character()),
      time_queued_sec = numeric(),
      time_running_sec = numeric(),
      worker_overhead_sec = numeric()
    ),
    R_runtime = 1L,
    n_cores = 1L,
    design = smoke_design
  )
  stopifnot(
    is.finite(runtime_summary$by_setting$mean_time_total_sec),
    isFALSE(runtime_summary$global$runtime_timestamps_complete)
  )
  message("Runtime smoke test passed.")
} else {
  message(
    "Model packages are unavailable. The runtime smoke test was skipped."
  )
}

# Synthetic setting-level times make the projection calculations transparent.
synthetic_by_setting <- tibble(
  setting_id = c(1L, 2L),
  mean_time_total_sec = c(2, 2)
)
repetitions <- 1000L
workers <- 8L
serial_seconds <- sum(
  synthetic_by_setting$mean_time_total_sec * repetitions
)
ideal_seconds <- serial_seconds / workers

complete_scheduler_summary <- tibble(
  runtime_timestamps_complete = TRUE,
  runtime_overhead_valid = TRUE,
  pilot_n_cores = workers,
  observed_overhead_factor = 1.25
)
complete_projection <- project_planned_study_runtime(
  synthetic_by_setting,
  complete_scheduler_summary,
  repetitions_per_setting = repetitions,
  n_cores = workers
)

stopifnot(
  abs(complete_projection$projected_serial_hours - serial_seconds / 3600) < 1e-12,
  abs(complete_projection$projected_ideal_parallel_hours - ideal_seconds / 3600) < 1e-12,
  abs(
    complete_projection$projected_adjusted_hours -
      (ideal_seconds * 1.25 / 3600)
  ) < 1e-12,
  isTRUE(complete_projection$runtime_timestamps_complete)
)

# Incomplete scheduler timestamps must not introduce an overhead factor into
# the projection. The internal-time projection remains available.
incomplete_scheduler_summary <- tibble(
  runtime_timestamps_complete = FALSE,
  observed_overhead_factor = 0.23
)
incomplete_projection <- project_planned_study_runtime(
  synthetic_by_setting,
  incomplete_scheduler_summary,
  repetitions_per_setting = repetitions,
  n_cores = workers
)

stopifnot(
  isFALSE(incomplete_projection$runtime_timestamps_complete),
  is.na(incomplete_projection$observed_overhead_factor_from_pilot),
  is.na(incomplete_projection$projected_adjusted_hours),
  abs(
    incomplete_projection$projected_ideal_parallel_hours -
      ideal_seconds / 3600
  ) < 1e-12
)



# Check internal task summaries and timestamp handling for a paused benchmark.
if (packages_available) {
  synthetic_results <- task_result[rep(1L, 4L), ]
  synthetic_results$task_id <- 1:4
  synthetic_results$run_id <- 1:4
  synthetic_results$time_total_sec <- rep(8, 4L)
  t0 <- as.POSIXct("2026-01-01", tz = "UTC")
  synthetic_jobs <- tibble(
    task_id = 1:4,
    submitted = t0 + c(0, 0, 0, 0),
    started = t0 + c(0, 0, 10, 10),
    done = t0 + c(10, 10, 20, 20),
    time_queued_sec = c(0, 0, 10, 10),
    time_running_sec = rep(10, 4L),
    worker_overhead_sec = rep(2, 4L)
  )
  execution <- list(
    schema = "runtime_execution_v1", submission_count = 1L,
    unknown_history = FALSE,
    profile = list(n_cores = 2L, cluster_function = "socket",
                   host = unname(Sys.info()[c("nodename", "sysname", "release", "machine")]))
  )
  continuous <- summarise_runtime_detail(
    synthetic_results, smoke_setting, synthetic_jobs, 4L, 2L,
    smoke_design, execution_metadata = execution
  )
  projected <- project_planned_study_runtime(
    continuous$by_setting, continuous$global, 100L, 2L
  )
  stopifnot(isTRUE(continuous$global$runtime_overhead_valid),
            abs(continuous$global$observed_overhead_factor - 1.25) < 1e-12,
            abs(projected$projected_adjusted_hours * 3600 - 500) < 1e-10)

  paused_jobs <- synthetic_jobs
  for (field in c("submitted", "started", "done")) {
    paused_jobs[[field]][3:4] <- paused_jobs[[field]][3:4] + 86400
  }
  resumed_execution <- execution
  resumed_execution$submission_count <- 2L
  paused <- summarise_runtime_detail(
    synthetic_results, smoke_setting, paused_jobs, 4L, 2L,
    smoke_design, execution_metadata = resumed_execution
  )
  paused_projection <- project_planned_study_runtime(
    paused$by_setting, paused$global, 100L, 2L
  )
  stopifnot(isTRUE(paused$global$runtime_timestamps_complete),
            isFALSE(paused$global$runtime_overhead_valid),
            is.na(paused$global$observed_overhead_factor),
            is.na(paused_projection$projected_adjusted_days),
            abs(paused_projection$projected_ideal_parallel_hours * 3600 - 400) < 1e-10,
            paused$global$pilot_makespan_sec == 86420)

  unverified <- summarise_runtime_detail(
    synthetic_results, smoke_setting, paused_jobs, 4L, 2L, smoke_design
  )
  stopifnot(isFALSE(unverified$global$runtime_overhead_valid),
            is.na(unverified$global$observed_overhead_factor))
  changed_workers <- project_planned_study_runtime(
    continuous$by_setting, continuous$global, 100L, 4L
  )
  stopifnot(is.na(changed_workers$projected_adjusted_days),
            abs(changed_workers$projected_ideal_parallel_hours * 3600 - 200) < 1e-10)
  partial <- summarise_runtime_detail(
    synthetic_results, smoke_setting, synthetic_jobs[1:3, ], 4L, 2L,
    smoke_design, execution_metadata = execution
  )
  stopifnot(isFALSE(partial$global$runtime_timestamps_complete))

  # Exercise stored history and the actual signature guard with a small
  # socket registry. Re-reading completed jobs does not add a submission;
  # retrying an unfinished job does, and changing workers is rejected.
  tiny_tasks <- make_tasks(smoke_setting, R = 2L, seed_base = 9000L)
  tiny_registry <- tempfile("runtime-history-")
  run_tiny <- function(workers = 1L) run_tasks(
    tiny_tasks, run_one_runtime_task, smoke_design,
    backend = "batchtools", registry_dir = tiny_registry,
    n_cores = workers, source_files = script_file,
    runtime_benchmark = TRUE
  )
  first <- run_tiny()
  first_summary <- summarise_runtime_detail(
    first$run_results, smoke_setting, extract_job_runtime_detail(first),
    2L, 1L, smoke_design, execution_metadata = first$runtime_execution
  )
  stopifnot(isTRUE(first_summary$global$runtime_timestamps_complete),
            isTRUE(first_summary$global$runtime_overhead_valid),
            is.finite(first_summary$global$observed_overhead_factor))
  repeated <- run_tiny()
  stopifnot(first$runtime_execution$submission_count == 1L,
            repeated$runtime_execution$submission_count == 1L,
            identical(first$run_results, repeated$run_results))
  mismatched <- try(run_tiny(2L), silent = TRUE)
  stopifnot(inherits(mismatched, "try-error"),
            grepl("execution_profile", as.character(mismatched), fixed = TRUE))
  batchtools::resetJobs(ids = data.table::data.table(job.id = 2L), reg = repeated$registry)
  resumed <- run_tiny()
  stopifnot(resumed$runtime_execution$submission_count == 2L)
  resumed_summary <- summarise_runtime_detail(
    resumed$run_results, smoke_setting, extract_job_runtime_detail(resumed),
    2L, 1L, smoke_design, execution_metadata = resumed$runtime_execution
  )
  stopifnot(isFALSE(resumed_summary$global$runtime_overhead_valid),
            is.na(resumed_summary$global$observed_overhead_factor))
  unlink(tiny_registry, recursive = TRUE)
}

message("All runtime-pilot tests passed.")
