# Check packages, protected files and the small test suites.
# No full pilot starts here.

source("R/workflow_helpers.R")
source("R/reproducibility.R")
check_reference_environment(packages = c(
  "batchtools", "dplyr", "tibble", "glmnet", "ranger", "snow", "ggplot2", "Matrix"
))
protected <- utils::read.csv("environment/protected_files.csv", stringsAsFactors = FALSE)
if (!all(file.exists(protected$path))) stop("A protected external file is missing.")
if (!identical(unname(tools::md5sum(protected$path)), protected$md5) ||
    any(file.info(protected$path)$size != protected$size_bytes)) {
  stop("An external file differs from its recorded bytes.")
}
scripts <- list.files(".", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
scripts <- scripts[!grepl("/(outputs|results)/", scripts)]
for (script in scripts) invisible(parse(file = script, keep.source = FALSE))
source("tests/run_all_tests.R")
cat("\nPROJECT VALIDATION COMPLETED.\n")
