# Full-grid counts, result selection and restart safeguards.
# The production pilots are not run by these tests.

source("R/planning_report.R")
source("R/elastic_net_intercept_check.R")
grid <- make_all_settings(make_design())
stopifnot(nrow(grid) == 192L, project_config$precision_R == 500L,
          project_config$rf_R == 100L,
          project_config$main_R == 5300L,
          identical(elastic_net_intercept_model_spec()$alpha, c(0.95, 1 / 3)))

precision <- list(
  settings = grid, design = make_design(),
  setting_summary = data.frame(setting_id = grid$setting_id, n_success = 500L,
                               R_required_sensitivity = c(5242L, rep(4300L, 191L))),
  global_summary = data.frame(n_tasks = 96000L, n_success = 96000L, n_failed = 0L,
                             max_sensitivity_tau = sqrt(5242) * 0.002)
)
check_precision_run(precision)
old <- precision
old$setting_summary$n_success <- 100L
stopifnot(inherits(try(check_precision_run(old), silent = TRUE), "try-error"))
subset <- precision
subset$settings <- head(grid, 30L)
stopifnot(inherits(try(check_precision_run(subset), silent = TRUE), "try-error"))
bad <- precision
bad$setting_summary$setting_id[2L] <- bad$setting_summary$setting_id[1L]
stopifnot(inherits(try(check_precision_run(bad), silent = TRUE), "try-error"))

candidates <- make_planning_candidates(precision)
stopifnot(!any(candidates$is_selected_main_R),
          candidates$repetitions_per_setting[candidates$scenario == "round_up_100"] == 5300L,
          candidates$repetitions_per_setting[candidates$scenario == "round_up_500"] == 5500L,
          !candidates$meets_pilot_point_target[candidates$scenario == "comparison_5200"])
selected <- make_planning_candidates(precision, main_R = project_config$main_R)
stopifnot(sum(selected$is_selected_main_R) == 1L,
          tail(selected$repetitions_per_setting, 1L) == 5300L,
          inherits(try(make_planning_candidates(precision, 1.5), silent = TRUE), "try-error"))

# Corrected families reserve separate beta, data, fold and RF seeds.
families <- c(1000, 300000000, 400000000, 500000000, 600000000, 650000000,
              1300000000, 1500000000, 1600000000, 1700000000)
lower <- families + 100000 + 11
upper <- families + 192 * 100000 + 9999 * 10 + 4
stopifnot(all(upper < .Machine$integer.max),
          all(head(upper, -1L) < tail(lower, -1L)))
tasks <- make_tasks(grid, R = 500L, seed_base = project_config$precision_seed_base)
stopifnot(nrow(tasks) == 96000L,
          !anyDuplicated(as.vector(outer(tasks$seed, 1:4, "+"))))

rf <- list(settings = grid, R_sensitivity = 100L,
           tree_counts = c(200L, 500L, 1000L),
           design = list(rf_tree_sensitivity_seed_policy = "separated"),
           setting_summary = expand.grid(setting_id = grid$setting_id,
                                         num_trees = c(200L, 500L, 1000L)))
rf$setting_summary$n_success <- 100L
check_rf_run(rf)
rf$design$rf_tree_sensitivity_seed_policy <- "old"
stopifnot(inherits(try(check_rf_run(rf), silent = TRUE), "try-error"))

# A completion record must match the sources, setup and saved output.
local({
  td <- tempfile("workflow-completion-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  source_file <- file.path(td, "source.R")
  writeLines("x <- 1", source_file)
  result_file <- file.path(td, "summary.rds")
  result <- list(value = 17L, saved_files = c(rds = result_file))
  saveRDS(result, result_file)
  spec <- list(R = 500L)
  mark_run_completed(result, td, spec, source_file)
  stopifnot(identical(load_completed_run(td, spec, source_file)$value, 17L),
            inherits(try(load_completed_run(td, list(R = 100L), source_file),
                         silent = TRUE), "try-error"))
  writeLines("x <- 2", source_file)
  stopifnot(inherits(try(load_completed_run(td, spec, source_file), silent = TRUE),
                    "try-error"))
  writeLines("x <- 1", source_file)
  saveRDS(list(value = 18L), result_file)
  stopifnot(inherits(try(load_completed_run(td, spec, source_file), silent = TRUE),
                    "try-error"))
  unlink(result_file)
  stopifnot(is.null(load_completed_run(td, spec, source_file)))
})

# Check the stage order and stop after a failed stage using substitutes.
local({
  code <- parse("run_planning_pilots.R", keep.source = FALSE)
  env <- new.env(parent = globalenv())
  env$source <- function(...) invisible(NULL)
  env$seen <- character()
  env$sys.source <- function(file, envir) {
    env$seen <- c(env$seen, file)
    invisible(NULL)
  }
  eval(code[[1L]], envir = env)
  env$run_planning_pilots()
  expected <- c("validate.R", "run/02_run_precision_pilot.R",
                "run/04_run_rf_tree_sensitivity.R", "run/10_summarize_planning.R",
                "run/11_export_results.R")
  stopifnot(identical(env$seen, expected))
  env$seen <- character()
  env$sys.source <- function(file, envir) {
    env$seen <- c(env$seen, file)
    if (grepl("02_run", file)) stop("Controlled precision failure.")
    invisible(NULL)
  }
  stopifnot(inherits(try(env$run_planning_pilots(), silent = TRUE), "try-error"),
            identical(env$seen, expected[1:2]))
})

cat("Workflow safeguards passed.\n")
