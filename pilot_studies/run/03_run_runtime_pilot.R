# Repeat the timing benchmark. The study plan uses main_R from run_config.R.

source("R/workflow_helpers.R")
source("R/runtime_pilot.R")
runtime_pilot <- run_runtime_pilot_8cores(
  output_dir = "outputs/runtime_pilot", R_runtime = 50L,
  n_cores = project_config$n_cores,
  registry_dir = "outputs/runtime_pilot/registry_R50"
)
writeLines(capture.output(sessionInfo()), file.path("outputs/runtime_pilot",
           paste0("session_info_", runtime_pilot$run_stamp, ".txt")))
print(runtime_pilot$runtime_global)
