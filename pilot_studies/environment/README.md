# Environment

Reference versions remain R 4.5.0, glmnet 4.1-8, ranger 0.17.0, batchtools 0.9.17, dplyr 1.1.4, tibble 3.2.1 and Matrix 1.7-3. Snow and ggplot2 are also needed. The optional full workflow needs xtable for the supplied binary script.

Missing packages stop execution. Version differences are reported as warnings. The standard starter installs nothing. If needed and internet access is available, run `source("run/00_install_packages.R")` separately.

Model threads are limited to one inside each R worker. The configured worker count is eight by default. The actual session is saved with new results.

Protected-file hashes refer to the preserved correlation/binary-use-case source files and included correlation figures. Package-file hashes describe this distribution and exclude the manifest itself.
