# Code guide

For first-time execution, read [RUN_GUIDE.md](RUN_GUIDE.md). It explains the console commands, dependencies, output folders and restart procedure. This document describes the implementation and scientific setup.

## Shared setup

`R/pilot_core.R` contains the existing simulation, fitting and batchtools infrastructure. `R/reproducibility.R` contains the version checks and seed arithmetic. The simulation design for the main study is:

| Factor            | Values            |
|-------------------|-------------------|
| n                 | 200, 1,000        |
| p/n               | 0.2, 5            |
| rho               | 0, 0.6, 0.8, 0.95 |
| a                 | 0.25, 0.75        |
| k                 | 1, 5, 10          |
| SNR               | 0.5, 2            |
| Block size        | 10                |
| Test observations | 5,000             |

This gives 192 settings. Elastic Net (EN) with alpha = 0.95 and 1/3 is used where EN is fitted. The precision pilot fits LASSO and RF only, and not EN or Ridge.

Here, n is the training-sample size, p/n controls dimensionality, rho is within-block predictor correlation, a is the fraction of active blocks and k is the number of active variables per active block. Gaussian outcomes follow Y = X beta + error with noise variance 1. Coefficients are scaled such that Var(X beta) / Var(error) equals the requested population SNR.

The population intercept is zero. The glmnet fits retain an unpenalized intercept and standardize predictors. Lambda is selected by five-fold cross-validation using MSE and predictions use lambda.min. Methods fitted within one repetition share its training/test data and CV folds. RF uses 200 trees, mtry = floor(sqrt(p)), min.node.size = 5 and one model thread, with tree counts varied in the RF sensitivity check.

For full reproducibiltity, RNG seeds are set in each repetition, calculated as: Task seed = seed_base + 100,000 \* setting_id + 10 \* run_id. Seeds using offsets +1, +2, +3 and +4 are set before generating beta, data, CV folds and fitting the RF, respectively. This prevents the reuse of component seeds across adjacent repetitions. Distinct seeds do not establish mathematical independence of whole random streams.

| Run                  |     Seed base |
|----------------------|--------------:|
| Precision            |         1,000 |
| Runtime              |   300,000,000 |
| RF sensitivity       |   400,000,000 |
| Test-set exact arm   |   500,000,000 |
| Test-set runtime arm |   600,000,000 |
| Test-set evaluation  |   650,000,000 |
| SNR check            | 1,300,000,000 |
| LASSO check          | 1,500,000,000 |
| EN check             | 1,600,000,000 |
| Ridge check          | 1,700,000,000 |

The precision pilot retains the seed family used in the R100 run, so its first 100 tasks per setting use the same seeds. The new run uses a separate R500 registry and calculates all 500 tasks per setting. RF sensitivity has its own non-overlapping seed family.

## Launchers

| Script | Work |
|----|----|
| `start.R` | Locate the project and refresh completed planning results |
| `update_planning.R` | Summary and result export without model fits |
| `run_planning_pilots.R` | Validation, run or reuse precision/RF pilots, summary, export |
| `resume.R` | Continue the precision/RF workflow only |
| `status.R` | Inspect precision/RF completion records or registry status and failed job IDs |
| `run_all.R` | All scripts 01 to 11 in run/, after validation |
| `run/00_install_packages.R` | Optional installation of missing packages |
| `run/01_run_testset_size_pilot.R` | Exact variance and test-set timing arms |
| `run/02_run_precision_pilot.R` | All 192 settings with R_sd = 500 |
| `run/03_run_runtime_pilot.R` | Runtime benchmark, 50 repetitions per setting |
| `run/04_run_rf_tree_sensitivity.R` | All 192 settings with 100 repetitions and three tree counts |
| `run/05_run_signal_sanity.R` | Full-grid SNR, known-coefficient prediction and intercept/model check, 50 repetitions |
| `run/06_run_lasso_intercept_check.R` | LASSO follow-up on the 12 selected settings, 50 repetitions |
| `run/07_run_elastic_net_intercept_check.R` | Same 12 settings, two EN alphas, 50 repetitions |
| `run/08_run_ridge_intercept_check.R` | Same 12 settings, Ridge alpha = 0, 50 repetitions |
| `run/09_run_correlation_binary_checks.R` | Correlation diagnostics and binary-use-case prevalence and RF summaries |
| `run/10_summarize_planning.R` | R and runtime tables from completed pilots |
| `run/11_export_results.R` | Review archive without execution registries |

The LASSO/EN/Ridge follow-up settings are IDs 145 to 156: n = 1,000, p = 5,000, rho = 0, and all a/k/SNR combinations. The two EN fits use the same generated data and CV folds within each task.

The source files for the analysis of the global correlation measures and the sensitivity analysis with binary outcome retain their supplied calculations, settings and seeds independent from the remaining pilots. The binary generator uses a logistic model with the same family of block-correlated predictors. The Gaussian seed and EN configuration applies to the modules in `R`. The existing correlation PDFs are retained.

## Restart and completion

The existing batchtools runner records task/design/source signatures, locks the registry during a run and checks for active jobs before restarting. Unfinished tasks can be restarted in the same registry after the previous session has ended. Failed jobs require `resubmit_errors = TRUE`.

`resume.R` resumes only the precision and RF stages. Runtime and diagnostic registries are continued through launchers 03 and 05 to 08. The runtime launcher fixes `resubmit_errors = FALSE`, so the configuration switch does not resubmit its failed jobs. Launchers 01 and 09 do not expose task-level registry continuation. See `RUN_GUIDE.md` for the individual commands.

`R/workflow_helpers.R` adds completion records for the precision and RF stages. A record includes the specification, source hashes and saved-output hashes. Missing output files cause the stage to collect its completed registry again. Changed setup, code or result bytes stop reuse.

`R/planning_report.R` requires 500 successful repetitions for every setting of the precision pilot, validates the RF counts and checks that the runtime grid matches. `main_R` only affects the report. Changing it does not invalidate the completed fitting stages.

## Agreed plan and source history

The agreed configuration for the main study derived from the precision pilot is `main_R = 5300L`: maximum calculated integer requirement 5242, rounded upward to 100. `outputs/planning_summary/main_study_plan.csv` records this choice.

The output-saving argument in the calculation functions is named `save_outputs`. File and function names use the project/workflow terminology. Calculations, fit settings and existing scientific result files are unchanged. Completion records retain the source hashes of the run and record the verified naming-only update. No source changes are accepted for completed-result reuse.

The result export excludes registries. If saved files are missing, recover them from the original working folder or recompute in a separate folder.
