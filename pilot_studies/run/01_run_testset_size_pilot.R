# Repeat the corrected test-set-size pilot if needed.

source("R/workflow_helpers.R")
source("R/testset_size_pilot.R")
testset_pilot <- run_testset_size_pilot(output_dir = "outputs/testset_size_pilot")
writeLines(capture.output(sessionInfo()), file.path("outputs/testset_size_pilot",
           paste0("session_info_", testset_pilot$run_stamp, ".txt")))
print(testset_pilot)
