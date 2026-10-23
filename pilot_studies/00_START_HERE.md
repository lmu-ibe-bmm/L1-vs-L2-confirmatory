# RIDGE–LASSO COMPARISON: START IN THE R CONSOLE

The pilot results and checks are already available. The agreed number of simulation repetitions, R, for the main simulation is 5300 per setting. The 500 repetitions refer to the precision pilot.

To get started without familiarity with the code: read `RUN_GUIDE.md`. The guide explains prerequisites, all starter scripts, result files, and how to resume individual checks. Basic knowledge of R and statistics is assumed.

1.  Extract the entire ZIP file.

2.  Run in R:

    ``` r
    source(file.choose())
    ```

3.  Select the file `start.R` in the `ridge-lasso-comparison` folder.

`start.R` sets the working directory and updates planning tables and the results archive from the existing results. It does not start any model fits.

## ALTERNATIVE

Open the project `ridge-lasso-comparison.Rproj` or set the project folder as the working directory. Then:

``` r
source("update_planning.R")
```

## MOST IMPORTANT RESULT FILES

-   `outputs/planning_summary/main_study_plan.csv`\
    Agreed number of simulation repetitions R = 5300, planning MCSE, and ideal runtime with eight workers.
-   `outputs/planning_summary/planning_status.txt`\
    Pilot maximum and decision rule for number of simulation repetitions.
-   `BINARY_USE_CASE.md`, `CORRELATION_MEASURES.md`\
    Extension of the Gaussian use case to binary outcomes, and analysis of overall correlation measures for real-data analysis.
-   `RESULTS.md`\
    Interpretation of the pilots, RF sensitivity, SNR and model checks.

## ADDITIONAL SCRIPTS

``` r
source("validate.R")
source("status.R")
```

`validate.R` checks code and small test cases. It does not start a full pilot. `status.R` shows existing completed results or the registry status.

Start the precision and RF pilots or use existing completed runs:

``` r
source("run_planning_pilots.R")
```

Run all planning pilots and checks:

``` r
source("run_all.R")
```

This is not required for the current analysis. It additionally starts test-set, runtime, SNR and model checks as well as correlation and binary checks. The confirmatory main simulation is not yet part of these scripts.

## RESUME THE PRECISION OR RF PILOT

Keep the entire working directory, including the registries. In the new R session, run:

``` r
source("status.R")
source("resume.R")
```

`resume.R` and `status.R` apply only to the precision and RF pilots, which are the most time-consuming. Completed tasks from these two pilots are reused. Do not start a second run in parallel in the same registry folder. Check failed tasks first. To resubmit them, set `resubmit_errors = TRUE` in `run_config.R`.

Instead, resume runtime, SNR and model checks using their respective starter scripts in `run/`. `RUN_GUIDE.md` names the exact files and explains the limitations of restarting.

## CONFIGURATION

`main_R = 5300`, `precision_R = 500`, `rf_R = 100` and `n_cores = 8`. If necessary, adjust `n_cores` to the workers actually available. Keep the documented pilot sizes and seed families. Elastic net uses `alpha = 0.95` and `1/3`. Seed comparisons are not included.
