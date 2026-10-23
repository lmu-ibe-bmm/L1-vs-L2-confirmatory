# Binary use case: extension to binary outcomes

This use case extends the study framework from prediction of a continuous Gaussian outcome to classification of a binary outcome. Both cases use the same family of block-correlated Gaussian predictors. The extension changes the conditional outcome distribution, the link between the linear predictor and the outcome, and the prediction-error measure.

The included binary analyses summarize outcome prevalence and RF classification performance for the binary extension. The planning pilots and confirmatory Ridge–LASSO design described elsewhere in this project concern the Gaussian outcome.

## From Gaussian regression to a binary outcome

In the Gaussian use case, Y = X beta + error, with independent Gaussian error and noise variance 1. In the binary use case:

-   Y given X follows a Bernoulli distribution with event probability pi.
-   logit(pi) = beta0 + X beta, with population intercept beta0 = 0.
-   The logistic transformation converts the linear predictor into a probability between zero and one.

The predictor remains linear in the coefficients, while the response is binary. Prediction is assessed here through RF out-of-bag classification error. The supplied script summarizes outcome prevalence and RF performance across the specified settings.

## Source files and execution

| Source file | Role |
|----|----|
| `R/run_logreg_pilots.R` | Generate binary outcomes and summarize prevalence and RF classification performance |

Run prevalence and RF prediction error checks together with correlation checks for real-world data analysis from the project directory with:

``` r
source("run/09_run_correlation_binary_checks.R")
```

The full launcher `run_all.R` includes this step. It requires `glmnet`, `ranger` and `xtable`. Binary summary tables are captured in `outputs/binary_check/binary_tables.tex`, with session information in the same folder.

## Binary analyses and output tables

The prevalence check records the proportion of observations with y = 1 in each generated dataset. The RF check records the out-of-bag classification error for each fitted forest (using default parameters of the ranger-implementation). Each LaTeX table reports the minimum, first quartile, median, mean, third quartile and maximum across repetitions, grouped by k and rho.

Both checks cover k = 1, 5 and 10 combined with rho = 0, 0.6, 0.8 and 0.95: 12 settings with n = 200, p = 1000, a = 0.25, b = 10 and SNR parameter 2. The prevalence check generates 1000 datasets per setting. The classification check fits only 100 forests per setting to keep the runtime low.
