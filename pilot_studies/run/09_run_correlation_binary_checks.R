# Run the scripts 

source("R/reproducibility.R")
check_reference_environment(packages = c("glmnet", "xtable", "ranger"))
dir.create("outputs/cor_summaries", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/binary_check", recursive = TRUE, showWarnings = FALSE)
message("Start recreating figures of correlation summaries.")
source("R/reproduce_figures_cor_summaries.R")
message("Finished recreating figures of correlation summaries.")
message("Start checks of binary DGM.")
utils::capture.output(
  source("R/run_logreg_pilots.R"),
  file = "outputs/binary_check/binary_tables.tex", split = TRUE
)
writeLines(capture.output(sessionInfo()), "outputs/binary_check/session_info.txt")
writeLines(capture.output(sessionInfo()), "outputs/cor_summaries/session_info.txt")
message("Finished checks of binary DGM.")
