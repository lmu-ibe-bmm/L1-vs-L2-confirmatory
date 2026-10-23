# Result guide

| Location | Contents |
|------------------------------------|------------------------------------|
| `outputs/precision_pilot_R500/` | Complete 192-setting precision pilot, 500 repetitions per setting |
| `outputs/rf_tree_sensitivity_R100/` | Corrected RF check, 100 repetitions and three tree counts |
| `outputs/planning_summary/main_study_plan.csv` | Selected R = 5300 and its planning/runtime estimates |
| `outputs/planning_summary/precision_runtime_candidates.csv` | Pilot requirement, rounding options and selected plan |
| `outputs/planning_summary/precision_global_summary.csv` | Counts, maximum tau and global requirements |
| `outputs/planning_summary/precision_top_settings.csv` | Settings with the largest estimated requirements |
| `outputs/planning_summary/rf_tree_count_summary.csv` | Paired tree-count sensitivity summaries |
| `results/` | Previously completed pilots and SNR/model checks |

The selected main-study number of repetitions, R, follows from the maximum number calculated in the complete R500 precision pilot rounded upward to 100. Candidate rows for other rounding choices are comparisons. They are not upper uncertainty bounds or alternative agreed decisions.

For details on the MCSE estimation, refer to the simulation protocol. The separate R100 RF check uses different simulated samples and is reported as sensitivity analysis to check the influence of the number of trees of the random forest on the MCSE estimation.

Runtime days mean 24 hours of continuous operation at the configured worker count. The included benchmark has 50 repetitions per setting. Its ideal projection for R = 5300 is about 14 days with eight workers. Shutdowns and additional server load increase elapsed calendar time.

Runtime projections use the recorded task times and worker count. Run the runtime benchmark with `source("run/03_run_runtime_pilot.R")`, followed by `source("update_planning.R")`.

The SNR check covers all 192 settings and evaluates whether the signal is strong enough that the performance of a random forest exceeds that of an intercept-only model. The LASSO, EN and Ridge follow-ups cover IDs 145 to 156 only, which are those scenarios for which the random forest performance did not exceed that of the intercept-only model.

## Reading the tables

| Field or quantity | Meaning |
|------------------------------------|------------------------------------|
| setting_id and design factors | Identify the scenario. A setting-summary row aggregates repetitions within that scenario |
| n_success / n_failed | Completed successful repetitions and failures. Precision requires 500 successes per setting, RF requires 100 and the diagnostic checks require 50 |
| sd_mse\_\* and cv_mse\_\* | Standard deviation of test MSE across repetitions and coefficient of variation (CV), i.e. SD divided by mean MSE. Similar CVs do not imply equal prediction errors |
| corr_lasso_rf_raw | Correlation of paired LASSO/RF test MSE values across repetitions |
| surrogate_tau_max | Estimated influence-value SD in the most demanding prespecified SD/effect scenario for that setting |
| R_required_sensitivity | Setting-specific integer requirement before rounding to 100, with the lower limit of 1000 |
| estimated_planning_mcse | Estimated tau / sqrt(R) under the pilot assumptions. The target is 0.002 |
| mean_analytic_snr | Population signal/noise ratio from the theoretical covariance. It should equal the setting's SNR |
| mean_empirical_snr_test | Signal/noise variance ratio observed in the finite test samples. It can fluctuate around the population value |
| mean\_\*\_gain_vs_intercept | Mean of repetition-level values 1 - model test MSE / baseline test MSE. For example 0.02 means a 2% relative improvement |
| \*\_gain_ci_lower / \*\_gain_ci_upper | Descriptive t-interval bounds for the mean repetition-level gain. A lower bound above zero indicates improvement over the training-mean baseline in that check |
| Fields containing oracle | Prediction using the known true coefficients, used as a signal reference without fitting a model |

The diagnostic gains are means of repetition-level ratios. The confirmatory Ridge–LASSO estimand uses a ratio of expected MSEs. These quantities serve different purposes. The model checks use separate seed families, so their separate gain summaries do not form a paired between-model contrast.

## Following a result back to its script

Launchers 01 to 09 write their outputs to the folders listed in RUN_GUIDE.md. SNR and model checks save task-level CSV files, setting summaries, overview tables, failure records and session information. The precision/RF public exports retain aggregated planning quantities and saved summary objects, while their individual execution caches are excluded.

Launcher 10 combines the completed precision/RF pilots with the latest compatible runtime benchmark into outputs/planning_summary/. Launcher 11 exports results/, outputs/ and run_config.R. update_planning.R calls 10 and 11. It updates the planning tables and archive, while each diagnostic check creates its own summaries.

The current R = 5300 projection is in outputs/planning_summary/main_study_plan.csv.

RESULTS.md interprets the included runs. It is a written report and does not change automatically when a check is rerun. After a new completed check, source("update_planning.R") refreshes the planning tables and includes the new output files in the archive.

The archive excludes batchtools registries and individual job caches. Keep the original working folder for task-level inspection or continuation of an incomplete run.

`BINARY_USE_CASE.md` explains the extension from a Gaussian response to a binary outcome, its prevalence/RF checks, the shared correlation figures and the coefficient construction and scaling.
