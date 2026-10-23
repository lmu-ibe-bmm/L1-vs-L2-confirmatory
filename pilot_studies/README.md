# Ridge–LASSO comparison

This project combines all pilot studies and analyses:

-   **Test set size pilot:** compare the finite sample variability of the prediction errors for different test set sizes and decide on the test set size to use for the main study and the remaining pilots

-   **Precision pilot**: calculate the required number of simulation repetitions to achieve the desired Monte Carlo Standard Error (MCSE) using a comparison of LASSO and a random forest

-   **Runtime pilot**: estimate the resulting overall runtime of the main study using a low number of simulation repetitions for each setting

-   **RF sensitivity analysis**: analyze how the number of trees affects the variability of the random forest and thus the results of the precision pilot

-   **SNR and model checks for the Ridge–LASSO study**: check whether the chosen signal-to-noise ratio levels are large enough that prediction models are able to achieve performance that is better than that of an intercept-only model

-   **Diagnostics of overall correlation measures**: calculate theoretical values of overall correlation measures used for the real-data analyses for the correlation levels of the planned study

-   **Checks preparing a binary-outcome use case**: prevalence and prediction errors of a random forest model for each setting in the planned small simulation study with binary outcome.

The completed results are included. The test set size is set to 5,000. The resulting number of repetitions for the main study is 5,300 repetitions per setting.

The complete precision pilot used 500 repetitions for each of the 192 settings. Its maximum integer requirement is 5,242. Rounding upward to the next 100 gives 5,300. The estimated maximum planning MCSE is 0.001989 against a target of 0.002.

## Start in R

Extract the whole ZIP. Open `ridge-lasso-comparison.Rproj`, or run `source(file.choose())` and select `start.R` inside the extracted folder. This locates the project and refreshes the planning tables and result archive without fitting models.

If the working directory already points to this folder:

``` r
source("update_planning.R")
```

The selected plan is in `outputs/planning_summary/main_study_plan.csv`. Start with [00_START_HERE.txt](00_START_HERE.txt). [01_RUN_GUIDE.md](01_RUN_GUIDE.md) explains prerequisites, the purpose and output of every launcher, completed-result reuse and restarting individual checks. It assumes basic R and statistics knowledge, without prior familiarity with this code. [07_RESULTS.md](07_RESULTS.md) explains the results and the choice of R.

## Project contents

| Folder or file | Purpose |
|------------------------------------|------------------------------------|
| `R/` | Shared simulation core, planning and diagnostic functions |
| `run/` | Numbered launchers for each pilot, check and export |
| `tests/` | Formula, fitting, seed, execution and restart checks |
| `external/` | Preserved source files for shared correlation diagnostics and the binary use case |
| `results/` | Previously completed corrected-code pilots and model checks |
| `outputs/` | Completed R500 precision pilot, corrected RF check and current plan |
| `environment/` | Reference environment, file inventory and validation records |
| `run_config.R` | Pilot sizes, worker count, seed families and main_R = 5300 |

## Reproduce or continue calculations

`source("run_planning_pilots.R")` validates the project, runs or reuses the complete precision and RF pilots, and updates the summary and export. The included completion records allow these two completed stages to be reused without new production fits.

`source("run_all.R")` also runs the test-set pilot, runtime benchmark, SNR and model checks, shared correlation diagnostics and binary-use-case checks. It starts substantially more computational work. The launchers in `run/` can be selected individually.

For an interrupted precision or RF pilot, use `source("status.R")` and `source("resume.R")` in the original server working folder. These two entry points cover only the precision/RF workflow. Continue the runtime, SNR and model checks through their individual launchers, as explained in `RUN_GUIDE.md`. Keep the working folder and its registries. The distributed result export does not include per-job caches or registries. No second run should be started against an active registry.

## Setup and scope

The folder contains the code and results for the planning pilots and checks discussed so far. The main confirmatory simulation and its hypothesis analysis have not been run by these scripts.

The binary use case extends the Gaussian response model to a Bernoulli outcome with a logistic link. [05_BINARY_USE_CASE.md](05_BINARY_USE_CASE.md) explains the extension, its prevalence and RF summaries, the shared correlation figures and the coefficient construction and scaling.
