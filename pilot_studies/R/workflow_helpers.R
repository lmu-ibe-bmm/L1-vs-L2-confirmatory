# Shared configuration checks and records for completed precision/RF stages.

source("run_config.R")

# Validate the shared configuration and limit BLAS threads inside each R process.
# Inputs: none, reads project_config. Returns invisible TRUE or stops.
setup_project <- function() {
  required <- c("run_config.R", "R/pilot_core.R", "R/precision_pilot.R")
  if (!all(file.exists(required))) stop("Set the working directory to ridge-lasso-comparison.")
  positive <- c("precision_R", "rf_R", "n_cores")
  for (name in positive) {
    x <- project_config[[name]]
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) {
      stop("Invalid run_config.R value: ", name)
    }
  }
  if (project_config$precision_R != 500L || project_config$rf_R != 100L) {
    stop("These launchers use 500 precision and 100 RF repetitions.")
  }
  if (!identical(project_config$precision_seed_base, 1000L) ||
      !identical(project_config$rf_seed_base, 400000000L)) {
    stop("Keep the recorded precision and RF seed families in run_config.R.")
  }
  if (!is.logical(project_config$resubmit_errors) ||
      length(project_config$resubmit_errors) != 1L ||
      is.na(project_config$resubmit_errors)) stop("Invalid resubmit_errors setting.")
  x <- project_config$main_R
  if (length(x) != 1L ||
      (!is.na(x) && (!is.numeric(x) || !is.finite(x) || x < 1 || x != floor(x)))) {
    stop("main_R must be NA or one positive integer.")
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
             NUMEXPR_NUM_THREADS = "1")
  return(invisible(TRUE))
}

# Fingerprint existing source files for completed-result reuse.
# Inputs: files. Returns an MD5 vector named by the supplied paths.
source_md5 <- function(files) {
  if (!all(file.exists(files))) stop("A required source file is missing.")
  return(setNames(unname(tools::md5sum(files)), files))
}

# Load a completed stage only when its specification, sources and saved files match.
# Inputs: output_dir, specification and source_files. Returns the saved summary.
# Returns NULL for a missing marker/output. Mismatched code or results cause an error.
load_completed_run <- function(output_dir, specification, source_files) {
  file <- file.path(output_dir, "run_complete.rds")
  if (!file.exists(file)) return(NULL)
  record <- readRDS(file)
  if (!identical(record$specification, specification) ||
      !identical(record$source_md5, source_md5(source_files))) {
    stop("The completed run does not match this setup: ", output_dir,
         ". Keep the existing folder and use a separate run folder.")
  }
  if (!all(file.exists(record$files))) return(NULL)
  if (!identical(unname(tools::md5sum(record$files)), record$file_md5)) {
    stop("A completed result was changed: ", output_dir)
  }
  message("Using the completed run in ", output_dir)
  return(readRDS(record$summary_file))
}

# Record the exact setup, sources and output hashes after a successful stage.
# Inputs: result, output_dir, specification and source_files.
# Writes run_complete.rds and returns its record invisibly.
mark_run_completed <- function(result, output_dir, specification, source_files) {
  files <- unname(result$saved_files)
  if (!length(files) || !all(file.exists(files))) stop("Result files are incomplete.")
  record <- list(
    specification = specification,
    source_md5 = source_md5(source_files),
    files = files,
    file_md5 = unname(tools::md5sum(files)),
    summary_file = unname(result$saved_files[["rds"]]),
    completed_at = Sys.time()
  )
  saveRDS(record, file.path(output_dir, "run_complete.rds"))
  return(invisible(record))
}

# Attach the actual software session to a newly completed result.
# Inputs: result and output_dir. Returns the result with its added session-file path.
save_run_session <- function(result, output_dir) {
  path <- file.path(output_dir, paste0("session_info_", result$run_stamp, ".txt"))
  writeLines(capture.output(sessionInfo()), path)
  result$saved_files <- c(result$saved_files, session_info = path)
  saveRDS(result, result$saved_files[["rds"]])
  return(result)
}

# Reject partial runs, the earlier R100 pilot and selected-setting subsets.
# Inputs: result and expected counts. Returns invisible TRUE for the complete R500 grid.
check_precision_run <- function(result, R = 500L, n_settings = 192L) {
  counts <- result$setting_summary$n_success
  g <- result$global_summary
  if (nrow(result$settings) != n_settings || nrow(result$setting_summary) != n_settings ||
      length(counts) != n_settings || anyNA(counts) || any(counts != R) ||
      g$n_tasks != n_settings * R || g$n_success != n_settings * R || g$n_failed != 0L ||
      anyDuplicated(result$setting_summary$setting_id) ||
      !setequal(result$settings$setting_id, result$setting_summary$setting_id) ||
      any(!is.finite(result$setting_summary$R_required_sensitivity))) {
    stop("The complete 192-setting precision pilot with 500 repetitions is required.")
  }
  return(invisible(TRUE))
}

# Require all settings and the three reviewed tree-count conditions.
# Inputs: result and expected counts. Returns invisible TRUE for the corrected R100 RF grid.
check_rf_run <- function(result, R = 100L, n_settings = 192L) {
  s <- result$setting_summary
  if (nrow(result$settings) != n_settings || result$R_sensitivity != R ||
      !identical(as.integer(result$tree_counts), c(200L, 500L, 1000L)) ||
      !identical(result$design$rf_tree_sensitivity_seed_policy, "separated") ||
      nrow(s) != n_settings * 3L || anyNA(s$n_success) || any(s$n_success != R) ||
      anyDuplicated(paste(s$setting_id, s$num_trees)) ||
      !setequal(s$setting_id, result$settings$setting_id) ||
      !setequal(s$num_trees, result$tree_counts)) {
    stop("The corrected full-grid RF tree-count run is incomplete.")
  }
  return(invisible(TRUE))
}

setup_project()
