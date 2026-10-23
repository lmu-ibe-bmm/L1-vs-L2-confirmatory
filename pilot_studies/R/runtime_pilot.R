###############################################################
# EIGHT-WORKER RUNTIME PILOT WITH DETAILED TIME COMPONENTS
###############################################################
# PURPOSE
# This script measures the computing time of one complete simulation
# repetition and supplies the measurements needed to project the resources for
# the planned study. It does not select the repetition count and its numerical
# model comparisons are not interpreted as scientific results.
#
# CONTRIBUTION TO THE OVERALL STUDY
# The precision pilot determines how many repetitions are needed, but it does
# not reveal whether that workload is operationally feasible. This runtime
# pilot measures the same full computational path that a study repetition uses
# and translates the observed setting-specific times into an eight-worker time
# projection. Keeping runtime planning separate prevents speed measurements
# from being mistaken for evidence about method performance.
#
# WHAT IS TIMED
# Every task generates training and independent test data, creates one shared
# set of cross-validation folds, and fits Ridge, LASSO, Elastic Net(1/3), and
# Elastic Net(0.95). Prediction, prediction mean squared error, prediction mean
# absolute error, and overall, active, and inactive coefficient mean squared
# errors are included in the measured time. Scheduler timestamps are recorded
# separately when batchtools provides them.
#
# METHODOLOGICAL BASIS
# Morris et al. (2019) and Siepe et al. (2024) support the explicit separation
# and reporting of simulation aims, data-generating mechanisms, methods, and
# performance measures. They do not prescribe this timing decomposition, the
# use of eight workers, or the scheduler-overhead adjustment. Those are
# study-specific planning choices reported here so that the feasibility
# calculation can be audited.
#
# RUNTIME PROJECTION
# For each setting, the pilot estimates the mean time of one complete task.
# The projected serial time is 5,200 times the sum of those setting-specific
# means. Dividing by eight gives the ideal parallel duration. When every
# scheduler timestamp is available and the history records one submission
# with a matching worker configuration, the ideal duration is multiplied by
# the observed pilot overhead factor. A resumed benchmark has no adjusted
# projection; its internal task-time projection remains available. The included
# reference benchmark produces approximately 12.19 adjusted days.
#
# One projected day is 24 hours of continuous operation at the stated worker
# count. Scheduled shutdowns and other periods without computation add calendar
# time. The projection does not specify a calendar completion date.
#
# Runtime values are hardware- and workload-specific. Another machine should
# reproduce the calculation and file structure, but it is not expected to
# reproduce the exact included timing values.
#
# WHEN AND HOW TO RUN THIS SCRIPT
# The file defines functions and does not start a pilot when sourced. From the
# top-level package directory, first run run/00_install_packages.R, then
# tests/run_all_tests.R, and finally run/03_run_runtime_pilot.R. The launcher
# uses eight workers and 50 repetitions per setting, matching the documented
# pilot configuration. To examine another server configuration, call
# run_runtime_pilot_8cores() explicitly with different n_cores or R_runtime.
# Results are written to outputs/runtime_pilot/ by the user-facing launcher.
# Every numerical result, projection, and session record from one execution
# uses the same timestamp, so repeated benchmarks remain distinguishable.
###############################################################

# PURPOSE
# Locate another R file stored beside this runtime module.
#
# WHY THIS FUNCTION EXISTS
# The module may be sourced by a launcher, run with Rscript, or opened
# interactively. These entry points need not share a working directory.
# Searching the actual script location first lets this file source the common
# core without relying on the package's parent-directory name.
#
# INPUT
# filename: File name to locate, normally "pilot_core.R".
#
# RETURN VALUE
# A normalized absolute path to the first existing candidate. The function
# stops with an informative error if no candidate exists.
find_sibling_file <- function(filename) {
  # Candidate 1: the directory of the script supplied through Rscript --file=.
  candidates <- character(0)

  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", cmd, fixed = TRUE)
  if (length(hit) > 0L) {
    this_file <- sub("--file=", "", cmd[hit[1]], fixed = TRUE)
    this_file <- normalizePath(this_file, winslash = "/", mustWork = FALSE)
    candidates <- c(candidates, file.path(dirname(this_file), filename))
  }

  # Candidate 2: the directory stored by source() in the active call stack.
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

  # Candidate 3: the current working directory, useful in interactive sessions.
  candidates <- unique(c(candidates, filename))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    stop(
      "File '", filename,
      "' was found neither beside the running script nor in the working directory."
    )
  }

  # Return an absolute path so later working-directory changes cannot alter it.
  return(normalizePath(existing[1], winslash = "/", mustWork = TRUE))
}

# Load the shared design, data-generating mechanism, model-fitting, timing, and
# execution functions. The core fixes basic linear algebra subprogram and
# model-internal thread counts before any computation, so those identical
# settings are not repeated in this module.
script_file <- find_sibling_file("pilot_core.R")
source(script_file)

# PURPOSE
# Execute the complete runtime pilot on the intended eight-worker setup and
# save both task-level measurements and setting-level summaries.
#
# WHY THIS FUNCTION EXISTS
# Computing time depends on the full repetition—not only on a single model
# fit. This runner therefore times data generation, four glmnet fits,
# prediction, and metric calculation under the same design as the later study.
# The function also retains scheduler timestamps, where available, so queuing
# and worker overhead can be distinguished from internal statistical work.
#
# INPUTS
# output_dir: Directory for the technical R data file and three
#   comma-separated values tables.
# R_runtime: Number of timing repetitions performed in every setting. The
#   documented pilot uses 50.
# n_cores: Number of socket workers. The documented pilot uses eight.
# registry_dir: Persistent batchtools registry. Reusing it resumes finished
#   tasks instead of recomputing them.
#
# RETURN VALUE
# A list containing the shared run_stamp, generated-file map, design, settings,
# task plan, task-level timings, batchtools status, persisted execution history,
# scheduler timing details,
# and both runtime summaries. The same object is saved as a timestamped R data
# file. PMSE and PMAE are calculated inside each timed fit but are not copied
# into this object. The local registry path is also omitted because the
# registry is execution state rather than a deposited result.
run_runtime_pilot_8cores <- function(
    output_dir = file.path("outputs", "runtime_pilot"),
    R_runtime = 50L,
    n_cores = 8L,
    registry_dir = file.path(output_dir, "registry")) {

  # Step 1: create the common factorial design and one deterministic task for
  # every setting-by-repetition combination.
  settings <- make_all_settings(design)
  tasks <- make_tasks(
    settings,
    R_runtime,
    # Keep timing tasks separate from precision and diagnostic seed families.
    seed_base = 300000000L
  )

  # Step 2: execute the plan through a resumable socket-worker registry. Every
  # task calls run_one_runtime_task(), which measures the full repetition.
  task_run <- run_tasks(
    tasks = tasks,
    task_function = run_one_runtime_task,
    design = design,
    backend = "batchtools",
    registry_dir = registry_dir,
    n_cores = as.integer(n_cores),
    resume = TRUE,
    overwrite_registry = FALSE,
    resubmit_errors = FALSE,
    wait = TRUE,
    source_files = script_file,
    runtime_benchmark = TRUE
  )

  # Step 3: separate scheduler timestamps from internal task timings, then
  # summarize task measurements and scheduler timestamps separately.
  job_detail <- extract_job_runtime_detail(task_run)
  summaries <- summarise_runtime_detail(
    run_results = task_run$run_results,
    settings = settings,
    job_detail = job_detail,
    R_runtime = as.integer(R_runtime),
    n_cores = as.integer(n_cores),
    design = design,
    execution_metadata = task_run$runtime_execution
  )

  # Step 4: combine the technical quantities required to audit or reanalyse
  # the benchmark. run_one_runtime_task() deliberately returns timings and
  # tuning information, not method-specific PMSE or PMAE values. The generated
  # support and block-sign lists are deterministic from the stored setting and
  # seed, and they are not needed for a runtime audit. Remove these list columns
  # from the public technical object while retaining the scalar diagnostics.
  public_run_results <- task_run$run_results %>%
    select(-active_blocks, -active_variables, -block_signs)

  output <- list(
    run_stamp = NA_character_,
    design = design,
    settings = settings,
    tasks = tasks,
    run_results = public_run_results,
    batchtools_status = task_run$status,
    runtime_execution = task_run$runtime_execution,
    job_runtime_detail = job_detail,
    runtime_by_setting = summaries$by_setting,
    runtime_global = summaries$global
  )

  # Step 5: save the technical object plus readable comma-separated values
  # extracts. A timestamp prevents one execution from overwriting another.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- make_timestamp()

  f_rds <- file.path(output_dir, paste0("runtime_pilot_", stamp, ".rds"))
  f_by_setting <- file.path(
    output_dir,
    paste0("runtime_by_setting_", stamp, ".csv")
  )
  f_global <- file.path(
    output_dir,
    paste0("runtime_global_", stamp, ".csv")
  )
  f_job_detail <- file.path(
    output_dir,
    paste0("job_runtime_detail_", stamp, ".csv")
  )

  # Store the shared timestamp and paths in the returned and saved object. The
  # launcher reuses run_stamp for the projection and software-session record.
  output$run_stamp <- stamp
  output$saved_files <- c(
    rds = f_rds,
    by_setting = f_by_setting,
    global = f_global,
    job_runtime_detail = f_job_detail
  )

  # Write the readable extracts first. Save the self-describing R data file only
  # after every referenced numerical table has been created successfully.
  write.csv(
    summaries$by_setting,
    f_by_setting,
    row.names = FALSE
  )
  write.csv(
    summaries$global,
    f_global,
    row.names = FALSE
  )
  write.csv(
    job_detail,
    f_job_detail,
    row.names = FALSE
  )
  saveRDS(output, f_rds)

  # Step 6: print the global summary and the 20 slowest settings. This gives the
  # person monitoring the server run an immediate diagnostic.
  print(summaries$global)
  print(
    summaries$by_setting %>%
      arrange(desc(mean_time_total_sec)) %>%
      select(
        setting_id, n, p, rho, a, k, snr, n_test,
        mean_time_total_sec, mean_time_ridge_fit_sec,
        mean_time_lasso_fit_sec, mean_time_en_fit_sec,
        mean_time_en095_fit_sec
      ) %>%
      head(20)
  )

  return(output)
}


###############################################################
# EXECUTION
###############################################################
# Sourcing this file intentionally does not start the runtime pilot. The
# recommended entry point is run/03_run_runtime_pilot.R. For a custom run,
# execute the function directly after package installation and tests:
#
# runtime_output <- run_runtime_pilot_8cores(
#   output_dir = file.path("outputs", "runtime_pilot"),
#   registry_dir = file.path("outputs", "runtime_pilot", "registry")
# )
###############################################################
