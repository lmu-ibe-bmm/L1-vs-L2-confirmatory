###############################################################
# TESTS FOR SOFTWARE-VERSION WARNINGS AND TASK SEED BOUNDARIES
###############################################################
# PURPOSE
# Verify the intended warning/error distinction and the numerical seed ranges.
# These tests use no model fits and start no pilot simulation.

source(file.path("R", "reproducibility.R"))

# A known mismatch remains a warning even inside the strict suite, where it is
# caught and inspected before R's warn=2 option could promote it to an error.
warning_capture <- new.env(parent = emptyenv())
warning_capture$condition <- NULL
comparison <- withCallingHandlers(
  check_reference_environment(
    packages = "stats", reference_versions = c(stats = "99.0.0"),
    reference_R = NULL
  ),
  reference_version_warning = function(condition) {
    warning_capture$condition <- condition
    invokeRestart("muffleWarning")
  }
)
observed_warning <- warning_capture$condition
stopifnot(
  inherits(observed_warning, "warning"),
  inherits(observed_warning, "reference_version_warning"),
  !comparison$matches_reference,
  comparison$installed == as.character(utils::packageVersion("stats"))
)

# Dotted and hyphenated patch versions represent the same package version.
matched <- check_reference_environment(
  packages = "Matrix",
  reference_versions = c(Matrix = gsub(
    "\\.([0-9]+)$", "-\\1", as.character(utils::packageVersion("Matrix"))
  )),
  reference_R = NULL
)
stopifnot(matched$matches_reference)
missing_package <- try(check_reference_environment(
  packages = "unavailable_simulation_test_package", reference_R = NULL
), silent = TRUE)
stopifnot(inherits(missing_package, "try-error"))

# Boundary repetitions in adjacent settings must reserve disjoint +1 to +4
# component seeds. Invalid ranges are rejected before conversion to integers.
boundary <- make_task_seeds(c(1L, 1L, 2L), c(1L, 9999L, 1L), 1000L)
components <- as.vector(outer(boundary, 1:4, "+"))
stopifnot(!anyDuplicated(components), all(is.finite(components)))
for (invalid_run in c(0, 1.5, 10000, 1e20, NA_real_)) {
  stopifnot(inherits(try(make_task_seeds(1L, invalid_run, 1000L),
                         silent = TRUE), "try-error"))
}
stopifnot(inherits(try(make_task_seeds(1L, 1L, .Machine$integer.max),
                       silent = TRUE), "try-error"))

# Seed families used by the planning pilots have disjoint complete setting
# ranges, including every supported repetition and component offset.
families <- c(
  precision = 1000, runtime = 300000000,
  testset_exact = 500000000, testset_runtime = 600000000,
  testset_evaluation = 650000000
)
lower <- families + 100000 + 10 + 1
upper <- families + 192 * 100000 + 9999 * 10 + 4
stopifnot(all(upper < .Machine$integer.max))
order_families <- order(lower)
stopifnot(all(upper[head(order_families, -1L)] <
  lower[tail(order_families, -1L)]))

cat("Software-environment and seed-boundary tests passed.\n")
