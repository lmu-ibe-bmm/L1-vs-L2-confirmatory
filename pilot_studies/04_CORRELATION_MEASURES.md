# Overall measures of correlation

The purpose of this check is to assess which global collinearity measure reflects the intended correlation levels and supports comparisons between simulation scenarios and real datasets.

## Source files and execution

| Source file | Role |
|----|----|
| `R/reproduce_figures_cor_summaries.R` | Assess global measures of the shared predictor-correlation structure |

Recreate figures together with running the binary use case from the project directory with:

``` r
source("run/09_run_correlation_binary_checks.R")
```

The full launcher `run_all.R` includes this step. It requires `glmnet`, `ranger` and `xtable`. Correlation figures are written to `outputs/cor_summaries/` with session information in the same folder.

## Correlation diagnostics

The correlation figures show summary measures of the predictor covariance matrix. They plot its determinant and condition number against within-block correlation rho at different dimensions p.

The plots show the determinant and condition number for each dimension p across the specified rho values. The two figures are `det_cor_mat.pdf` and `cond_cor_mat.pdf`.
