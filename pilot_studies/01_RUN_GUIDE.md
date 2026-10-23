# Running Scripts: Getting Started and Workflow

This guide is intended for people with basic knowledge of R and statistics who are not yet familiar with the project code. It explains what the scripts calculate, how to start them, and where the results are located. The statistical design is described in [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md), and the analysis in [RESULTS.md](RESULTS.md).

## What Is Already Available

The folder contains the completed planning pilots and the corrected SNR and model checks. No new full pilot run is necessary for the current analysis.

The complete precision pilot yields a maximum integer requirement of 5,242. Rounding up to the next 100 gives 5,300. These scripts do not yet run the confirmatory main simulation or its hypothesis analysis.

## Prerequisites and First Start

1.  Extract the **entire ZIP file**. The subfolders must remain together.

2.  Use R. The production runs used R 4.5.0. The required packages are `glmnet`, `ranger`, `batchtools`, `snow`, `dplyr`, `tibble`, `ggplot2`, and `Matrix`. `xtable` is additionally required for the tables of the binary extension and `run_all.R`. The reference versions are listed in [environment/README.md](environment/README.md).

3.  Run the following in the R console:

    ``` r
    source(file.choose())
    ```

    In the file dialog, select `start.R` in the extracted project folder. The script sets the working directory and updates the planning tables and results export from the existing pilot results. **No models are refitted** in the process.

Alternatively, open `ridge-lasso-comparison.Rproj` in RStudio and then run `source("update_planning.R")`. Without RStudio, the working directory can be set to the extracted project folder using `setwd(...)`.

All commands below assume that `getwd()` points to this project folder. The computation starters do not install packages. Missing packages must be available beforehand. If internet access and a writable package library are available, missing packages can be installed separately:

``` r
source("run/00_install_packages.R")
```

On a server without internet access, the packages must already be installed. Missing packages cause the run to stop. Differences from the documented software versions are reported as warnings; they are not evidence of a faulty run. Identical results across different software versions are therefore not guaranteed.

## Which Command Fits Which Goal?

| Goal | Command | What happens? |
|------------------------|------------------------|------------------------|
| Update the existing plan | `source("update_planning.R")` | Recreates planning tables and the results archive without model fits |
| Check code and small test cases | `source("validate.R")` | Checks packages, protected files, syntax, and eleven test suites, including small real model fits and socket tests |
| Check completion or job status | `source("status.R")` | Checks existing precision/RF completions or their registries without starting jobs |
| Run the precision/RF workflow | `source("run_planning_pilots.R")` | Validates, runs 02 and 04 or reuses their valid completions, then runs 10 and 11 |
| Run all pilots and checks | `source("run_all.R")` | Validates and runs 01 through 11. Existing completions of 02/04 are reused; other stages may start new computations |
| Resume the precision/RF workflow | `source("resume.R")` | Starts the same workflow as `run_planning_pilots.R` in the existing working folder |

The included completion files allow the completed precision and RF pilots to be reused. Therefore, `run_planning_pilots.R` does not recompute these two pilots in the unchanged package. Validation nevertheless performs small test fits.

`run_all.R` is not needed to read the existing results and may initiate substantial computation. A completely new run of 02 or 04 requires a **separate working copy without their previous outputs and completion files**. For an entirely new run, the whole `outputs/` folder can be omitted from that separate copy. Keep the existing working folder and its results.

## All Numbered Starters

A setting is a combination of design factors; a repetition is a newly generated dataset with its associated model fits. The starters are located in `run/`. The files in `R/` define the functions and are not the entry points for complete computations.

Each starter is executed from the project folder using `source("run/DATEINAME.R")`. Example for the SNR check:

``` r
source("run/05_run_signal_sanity.R")
```

| Starter in `run/` | Purpose and scope | Results location |
|------------------------|------------------------|------------------------|
| `00_install_packages.R` | Installs only missing packages. Does not run a simulation | R package library |
| `01_run_testset_size_pilot.R` | Examines variability across different test set sizes and measures the computational cost of test evaluation | `outputs/testset_size_pilot/` |
| `02_run_precision_pilot.R` | LASSO/RF pilot across 192 settings × 500 repetitions, totaling 96,000 tasks. Estimates the required R for the main simulation | `outputs/precision_pilot_R500/` |
| `03_run_runtime_pilot.R` | Measures runtimes across 192 settings × 50 repetitions. Basis for extrapolation to the main simulation R | `outputs/runtime_pilot/` |
| `04_run_rf_tree_sensitivity.R` | Compares 200, 500, and 1,000 RF trees on the same data for each task, across 192 settings × 100 repetitions | `outputs/rf_tree_sensitivity_R100/` |
| `05_run_signal_sanity.R` | Checks the simulated signal, predictions using known coefficients, and RF against constant prediction using the training mean. 192 settings × 50 repetitions | `outputs/signal_sanity/` |
| `06_run_lasso_intercept_check.R` | Checks LASSO and intercept diagnostics on twelve selected settings, with 50 repetitions each | `outputs/lasso_intercept_check/` |
| `07_run_elastic_net_intercept_check.R` | Corresponding EN check with alpha = 0.95 and 1/3, jointly generated data and CV folds per task. Twelve settings × 50 repetitions | `outputs/elastic_net_intercept_check/` |
| `08_run_ridge_intercept_check.R` | Corresponding Ridge check with alpha = 0. Twelve settings × 50 repetitions | `outputs/ridge_intercept_check/` |
| `09_run_correlation_binary_checks.R` | Checks the joint correlation structure and the extension of the use case to binary outcomes through prevalence and RF classification checks. The existing code is retained | `outputs/cor_summaries/` and `outputs/binary_check/` |
| `10_summarize_planning.R` | Reads completed pilots and creates R, RF, and runtime tables. No new model fits | `outputs/planning_summary/` |
| `11_export_results.R` | Packages public results and configuration into an archive. No new model fits | `ridge-lasso-comparison-results.tar.gz` |

The repetition counts serve different purposes:

| Number | Meaning |
|------------------------------------|------------------------------------|
| 500 per setting | Precision pilot: estimates the number of repetitions required for the main simulation across all 192 settings |
| 100 per setting | RF sensitivity check with 200, 500, and 1,000 trees |
| 50 per setting | Runtime pilot and SNR/model checks, each within its documented scope |
| 5,300 per setting | Agreed number of repetitions for the subsequent main simulation |

The twelve settings in 06 through 08 are IDs 145 through 156: n = 1,000, p = 5,000, rho = 0, and all combinations of a, k, and SNR. These model checks therefore cover a selected subset. SNR check 05 covers all 192 settings.

Starters 01 through 09 compute their respective pilots, model checks, and use-case checks. They do not require all lower-numbered scripts to have been run beforehand. For 02/04, the reuse of existing completions described above applies.

Starter 10 requires the complete R500 precision pilot, the corrected R100 RF check, and a suitable runtime benchmark. A new benchmark from `outputs/runtime_pilot/` is preferred; otherwise, the included benchmark from `results/runtime_pilot/` is used. Starter 11 requires an existing planning summary. `update_planning.R` runs 10 and 11 in the correct order.

The joint correlation diagnostics and binary-extension checks are started via 09 so that their results are placed in the designated folders. [BINARY_USE_CASE.md](BINARY_USE_CASE.md) explains the transition from the Gaussian regression case to a Bernoulli outcome with a logistic link and contextualizes the checks and plots. The tests are started via `validate.R`. It is not necessary to run every file in `tests/` individually.

## Finding Results and Checking a Run

-   `results/` contains the already completed runs retained in the package. This serves as an archive of the results used as a basis for the simulation design to keep for comparison if the code is rerun.
-   `outputs/` at the moment also contains completed runs retained in the package. New runs also write here. Delete the contents before starting a rerun.
-   `outputs/planning_summary/main_study_plan.csv` contains the agreed R = 5,300 as well as the planning and runtime estimates.
-   `outputs/planning_summary/planning_status.txt` explains the pilot maximum and the decision rule.
-   The individual results or outputs folders contain CSV tables, RDS objects, and, depending on the check, figures and session information. [RESULTS_GUIDE.md](RESULTS_GUIDE.md) identifies the most important files.

CSV files can be opened in R or in a spreadsheet application. For example, the selected plan can be displayed with:

``` r
plan <- read.csv("outputs/planning_summary/main_study_plan.csv")
print(plan)
```

`RESULTS.md` is the written analysis of the included runs. It is not automatically rewritten by rerunning computations.

`RESULTS_GUIDE.md` explains the most important columns and their meanings. SNR and model checks generate their own summaries in their respective results folders. After a newly completed check, update the planning and results export using `source("update_planning.R")`. The new check analysis remains in its own tables.

For 02/04, the completion check and `status.R` confirm that the full results are available. When reusing a run, the console displays `Using the completed run in ...`. Other checks report completion in the console and save summary/error tables in their respective folders. An interrupted run or a table with fewer successful repetitions than planned does not yet constitute a complete analysis. Failed jobs must be examined first. The planning summary rejects incomplete precision/RF results.

## Resuming After an Interruption

Resuming requires the **original working folder with its registries**. These store completed tasks and the code/design state. The distributed ZIP and the results export do not contain individual job caches or registries.

After the previous R session ends, open a new session, set the working directory to the same project folder, and run the appropriate starter:

| Interrupted run | Command to resume |
|------------------------------------|------------------------------------|
| Precision or RF pilot | `source("status.R")`, then `source("resume.R")` |
| Runtime pilot | `source("run/03_run_runtime_pilot.R")` |
| SNR check | `source("run/05_run_signal_sanity.R")` |
| LASSO check | `source("run/06_run_lasso_intercept_check.R")` |
| EN check | `source("run/07_run_elastic_net_intercept_check.R")` |
| Ridge check | `source("run/08_run_ridge_intercept_check.R")` |
| Reproduce correlation plots and binary check | `source("run/09_run_correlation_binary_checks.R")` |

Registered, completed tasks are reused. `status.R` and `resume.R` cover only the precision and RF pilots, not all model checks.

Do not start a second computation in parallel on the same registry. Do not replace the code of a registry that has not yet completed with a different package version, because the signatures also check the source files.

Failed tasks are not resubmitted by default. After examination, `resubmit_errors = TRUE` can be set in `run_config.R` for 02, 04, and 05 through 08. Runtime starter 03 uses a fixed `resubmit_errors = FALSE`; this configuration change is therefore not sufficient for errors in the runtime pilot. First determine the specific cause of the error.

Starters 01 and 09 do not support resuming through a task registry. After an interruption, they are restarted and may repeat their computations.

## Configuration and Runtime

`run_config.R` specifies `main_R = 5300L`, `precision_R = 500L`, `rf_R = 100L`, and eight workers by default. Workers are R processes running in parallel. `n_cores` can be adjusted to the number of workers actually available. Pilot sizes and seed families are part of the specified setup and should be retained to reproduce it. Model threads are limited to one within each worker.

The current ideal runtime extrapolation for the main simulation is approximately 14 days with eight workers operating continuously, 24 hours a day. This is neither the duration of an individual pilot nor a guarantee of elapsed calendar time. Shutdowns and additional server load extend it.
