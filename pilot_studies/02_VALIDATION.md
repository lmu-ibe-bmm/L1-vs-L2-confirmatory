# Validation

Run `source("validate.R")` from the project folder. The planning and full-calculation starters call it before long runs. `start.R` and `update_planning.R` only refresh the included completed results.

Validation checks package availability, parses the R scripts, verifies the preserved correlation/binary-use-case source bytes and runs the shared-core, planning and diagnostic tests. The tests cover formulas, small real model fits, socket execution, registry safeguards, seed ranges, full-grid task counts, completed-result reuse and stage order. No complete pilot study run is started by validation.

Missing packages and unexpected warnings stop tests. Differences from the recorded R/package versions produce a named reference-version warning and allow execution. Session records retain the actual run environment.

The completed long pilots were run on an server with R 4.5.0.

The assembly audit verifies every source/launcher reference, all file and identifier names, computational equivalence, unchanged scientific result files, selected R = 5300, updated planning tables and package extraction.

`environment/package_files.csv` records the distribution's file sizes and MD5 hashes. `environment/protected_files.csv` is enforced for the distributed correlation/binary-use-case source files. Configuration changes for worker count, error resubmission and the agreed main-study R are supported.
