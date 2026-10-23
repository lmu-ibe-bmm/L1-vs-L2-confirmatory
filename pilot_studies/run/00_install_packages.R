# Install missing packages only. Existing versions are kept.

packages <- c("glmnet", "ranger", "batchtools", "snow", "dplyr",
              "tibble", "ggplot2", "Matrix", "xtable")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!length(missing)) cat("All required packages are already installed.\n")
