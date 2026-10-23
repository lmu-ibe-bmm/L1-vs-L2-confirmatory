###############################################################
# SHARED SOFTWARE-ENVIRONMENT AND SEED CHECKS
###############################################################
# PURPOSE
# Keep reference-version checks and task-seed arithmetic consistent across
# the shared core, the standalone test-set pilot, and diagnostic modules.
# Only base R is needed to load this file. Sourcing it starts no simulation.

# PURPOSE
# Return the software-version specification for the pilot environment.
#
# WHY IT EXISTS
# Reference versions are provenance information, rather than a requirement
# that prevents execution in a different compatible software environment.
#
# INPUTS
# None.
#
# RETURN VALUE
# A named character vector of recorded package versions.
reference_package_versions <- function() {
  return(c(
    glmnet = "4.1.8", ranger = "0.17.0", batchtools = "0.9.17",
    dplyr = "1.1.4", tibble = "3.2.1", Matrix = "1.7-3"
  ))
}

# PURPOSE
# Check availability and report differences from the reference environment.
#
# WHY IT EXISTS
# Missing packages prevent a run, whereas a different installed version should
# produce a warning and remain visible in the saved session information.
# A specific warning class lets the strict test runner retain this distinction.
#
# INPUTS
# packages: Packages needed by the current entry point.
# reference_versions: Named versions used for comparison.
# reference_R: Recorded R version, or NULL to omit the R comparison.
#
# RETURN VALUE
# A data frame of installed and reference versions, invisibly. Missing packages
# cause an error. Version differences produce a reference_version_warning.
check_reference_environment <- function(
    packages = names(reference_package_versions()),
    reference_versions = reference_package_versions(),
    reference_R = "4.5.0") {
  available <- vapply(
    packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
  )
  if (!all(available)) {
    stop("Missing R packages: ", paste(packages[!available], collapse = ", "))
  }

  installed <- vapply(packages, function(package) {
    return(as.character(utils::packageVersion(package)))
  }, character(1))
  reference <- unname(reference_versions[packages])
  comparison <- data.frame(
    package = packages, installed = unname(installed), reference = reference,
    stringsAsFactors = FALSE
  )
  if (!is.null(reference_R)) {
    comparison <- rbind(
      data.frame(
        package = "R", installed = as.character(getRversion()),
        reference = reference_R, stringsAsFactors = FALSE
      ), comparison
    )
  }

  comparison$matches_reference <- vapply(seq_len(nrow(comparison)), function(i) {
    if (is.na(comparison$reference[i])) return(NA)
    # package_version() treats dotted and hyphenated patch numbers equally.
    return(identical(
      package_version(comparison$installed[i]),
      package_version(comparison$reference[i])
    ))
  }, logical(1))
  different <- which(!is.na(comparison$matches_reference) &
    !comparison$matches_reference)
  if (length(different) > 0L) {
    details <- paste0(
      comparison$package[different], " installed ",
      comparison$installed[different], ", reference ",
      comparison$reference[different]
    )
    condition <- structure(
      list(
        message = paste(
          "Software versions differ from the recorded reference environment:",
          paste(details, collapse = ", "),
          "Execution continues. Record the current session information."
        ), call = NULL
      ),
      class = c("reference_version_warning", "warning", "condition")
    )
    warning(condition)
  }
  return(invisible(comparison))
}

# PURPOSE
# Generate the collision-free base seed for each setting-repetition pair.
#
# WHY IT EXISTS
# A gap of ten reserves separate component seeds +1 beta, +2 data, +3 CV folds,
# and +4 Random Forest. The setting gap permits at most 9,999 repetitions.
# Calculating with doubles before integer conversion prevents silent overflow.
# Distinct numeric seeds prevent this particular reuse, they do not prove
# independence of entire pseudorandom streams.
#
# INPUTS
# setting_id, run_id: Equally sized vectors of positive integer identifiers.
# seed_base: One non-negative integer identifying the seed family.
#
# RETURN VALUE
# Integer task seeds, or an error for invalid identifiers or an unsafe range.
make_task_seeds <- function(setting_id, run_id, seed_base) {
  if (length(seed_base) != 1L || !is.finite(seed_base) || seed_base < 0 ||
    seed_base != floor(seed_base)) {
    stop("seed_base must be one non-negative integer.")
  }
  if (length(setting_id) < 1L || length(setting_id) != length(run_id) ||
    any(!is.finite(setting_id)) || any(setting_id < 1) ||
    any(setting_id != floor(setting_id)) || any(!is.finite(run_id)) ||
    any(run_id < 1) || any(run_id > 9999) || any(run_id != floor(run_id))) {
    stop("Use positive integer setting IDs and run IDs from 1 to 9999.")
  }
  seed <- as.double(seed_base) + 100000 * as.double(setting_id) +
    10 * as.double(run_id)
  if (any(seed + 4 > .Machine$integer.max)) {
    stop("The task and component seeds exceed the R integer range.")
  }
  return(as.integer(seed))
}
