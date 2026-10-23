# Ridge–LASSO comparison: Results and choice of R

## 1. Agreed number of repetitions: R = 5300

For the main simulation, 5300 repetitions per setting are used. The decision follows the agreed planning rule: the maximum estimated requirement is determined across all 192 settings and the specified SD/effect scenarios of the complete precision pilot, then rounded up to the next 100.

The complete pilot with 500 repetitions each yields a maximum integer requirement of 5242. Thus ceiling(5242 / 100) × 100 = 5300. The maximum estimated influence-value SD is 0.14478995. With 5300 repetitions, the estimated maximum planning MCSE is 0.14478995 / sqrt(5300) = 0.00198884. It is therefore below the target value of 0.002.

The planning already accounts for SD multipliers up to 1.3 and assumed relative effects of -0.02, 0, and 0.02. The maximum occurs at multiplier 1.3 and effect -0.02. The choice of 5300 thus uses the maximum of the sensitivity grid.

5300 is a traceable point estimate for planning under the documented assumptions. The estimate varies between pilot samples. It does not guarantee that every subsequent pilot sample would yield a requirement below 5300. The separate RF check with 100 repetitions, for example, yields 5535 for 200 trees. This deviation is reported as remaining uncertainty in the planning estimate. The main-simulation R is determined from the complete 500-repetition precision pilot.

The planning values are based on RF/LASSO dependence and specified SD/effect scenarios. They are estimated planning MCSEs, not MCSEs already measured for the actual Ridge–LASSO contrast.

## 2. Completeness and implementation

The project contains all discussed Gaussian planning pilots and checks in one folder: test-set pilot, precision pilot, runtime pilot, RF sensitivity, SNR check, LASSO, EN, and Ridge checks. In addition, it includes the joint correlation diagnostics and the checks for extending the use case to binary outcomes. BINARY_USE_CASE.md describes this extension using a Bernoulli distribution with a logistic link.

| Component | Scope | Available results |
|------------------------|------------------------|------------------------|
| Complete precision pilot | 192 settings × 500 repetitions | 96000 successful, 0 failed |
| RF check | 192 settings × 100 repetitions, 3 RF variants each | 19200 tasks successful, 576 setting/tree-count rows |
| Runtime pilot | 192 settings × 50 repetitions | 9600 successful |
| Test-set pilot | Analytical precision arm and runtime arm | Results carried over |
| SNR check | 192 settings × 50 repetitions | 9600 successful |
| LASSO and Ridge checks | 12 settings × 50 repetitions each | 600 successful each |
| EN check | 12 settings × 50 repetitions × 2 alphas | 1200 successful model rows |

The project modules use fixed seeds. The task seed is seed_base + 100000 × setting_id + 10 × run_id. The components beta, data, CV folds, and RF use offsets +1, +2, +3, and +4. EN uses alpha = 0.95 and 1/3.

The calculations for the correlation diagnostics and binary extension remain unchanged. The binary data generator retains its previous coefficient, scaling, and seed settings. The uniform Gaussian configuration applies to the Gaussian project modules.

## 3. Precision pilot

The design comprises n = 200 and 1000, p/n = 0.2 and 5, rho = 0, 0.6, 0.8, and 0.95, a = 0.25 and 0.75, k = 1, 5, and 10, and SNR = 0.5 and 2. This results in 192 settings. Block size b = 10 and n_test = 5000 remain constant.

| Setting |   n |    p |  rho |    a |   k | SNR | R requirement |
|--------:|----:|-----:|-----:|-----:|----:|----:|--------------:|
|      50 | 200 | 1000 |    0 | 0.25 |   1 |   2 |          5242 |
|      62 | 200 | 1000 |  0.6 | 0.25 |   1 |   2 |          3736 |
|      96 | 200 | 1000 | 0.95 | 0.75 |  10 |   2 |          3650 |
|      94 | 200 | 1000 | 0.95 | 0.75 |   5 |   2 |          3423 |

The maximum is determined by setting 50. Higher requirements are concentrated at n = 200, p/n = 5, and SNR = 2. Neither high dimensionality nor high correlation alone determines the requirement.

| SD multiplier | Maximum requirement across settings and effect scenarios | Rounded up to 100 |
|-----------------------:|-----------------------:|-----------------------:|
| 1.0 | 3856 | 3900 |
| 1.1 | 4259 | 4300 |
| 1.2 | 4721 | 4800 |
| 1.3 | 5242 | 5300 |

## 4. RF sensitivity

Within each task, the RF check uses the same data and the same LASSO fit for the three tree counts. The check has its own seed range and 100 repetitions per setting.

| Trees | Median SD_RF / SD_LASSO | Median RF-LASSO correlation | Maximum R requirement | Rounded up to 100 |
|--------------:|--------------:|--------------:|--------------:|--------------:|
| 200 | 0.948 | 0.625 | 5535 | 5600 |
| 500 | 0.941 | 0.630 | 5261 | 5300 |
| 1000 | 0.938 | 0.629 | 5337 | 5400 |

Compared with 200 trees, the RF error SD decreases by around 0.47 percent in the median with 500 trees and by around 0.53 percent with 1000 trees. The median paired change in correlation is approximately +0.0041 and +0.0046, respectively. The differences are small.

More trees do not lead to a strictly monotonically decreasing R requirement in this finite pilot sample. The sensitivity check provides no clear reason to change the 200 trees of the auxiliary RF pilot. The RF tree count concerns the auxiliary model used for planning.

The separate 100-repetition check is smaller than the complete 500-repetition precision pilot and uses different simulated samples, resulting in different R requirements. The R requirements from the precision pilot are kept as they are calculated using a higher number of iterations and therefore present a less variable estimate.

## 5. SNR and model checks

The difference between the analytical SNR in the simulation and the target SNR is at most approximately 2.2 × 10\^-15. The data generator thus sets the intended population signal correctly to numerical precision. Signal strength in individual finite samples still varies.

In the complete RF-SNR check, 187 of 192 descriptive intervals for the mean relative MSE improvement over the constant prediction using the training mean lie entirely above zero. At SNR = 0.5, this applies to 91 of 96 settings; at SNR = 2, to all 96.

The median relative improvement is 6.50 percent at SNR = 0.5 and 14.89 percent at SNR = 2. The five settings without RF improvement have n = 1000, p = 5000, rho = 0, and SNR = 0.5. Their negative mean values are small, approximately -0.076 to -0.135 percent.

The in-depth LASSO, EN, and Ridge checks concern only the 12 settings with n = 1000, p = 5000, and rho = 0 across all a/k/SNR combinations.

| Model | Median improvement at SNR = 0.5 | Median improvement at SNR = 2 | Settings with interval entirely above zero |
|------------------|-----------------:|-----------------:|-----------------:|
| LASSO | -0.074 percent | 1.175 percent | 5 of 12 |
| EN alpha = 0.95 | -0.050 percent | 1.045 percent | 5 of 12 |
| EN alpha = 1/3 | -0.029 percent | 1.489 percent | 6 of 12 |
| Ridge | 1.865 percent | 8.093 percent | 12 of 12 |

Positive values mean lower test MSE than prediction using the training mean. Ridge shows consistently positive improvements in the 12 settings examined, and is thus able to pick up the signal as intended.

The signal is present and the models exploit it differently. The difficult settings remain in the design. These checks use their own seed ranges and are descriptive diagnostics of learnability. They do not constitute a paired Ridge–LASSO contrast and do not test a confirmatory hypothesis.

## 6. Test set and runtime

n_test = 5000 is chosen. In the Gaussian model used, the CV of the mean squared prediction errors with approximately centered errors is roughly sqrt(2 / n_test). This yields about 2 percent with 5000 test observations. Similar CVs across models are therefore explainable and do not imply identical prediction errors or an absence of signal.

The analytical values at n_test = 2500, 5000, 10000, and 20000 are approximately 2.83, 2.00, 1.41, and 1.00 percent. The empirical SDs in the separate runtime arm are based on only five repetitions and provide less stable estimates.

The runtime projection is based on the available corrected benchmark with 192 settings and 50 repetitions each. For R = 5300, eight workers yield an ideal parallel runtime of approximately **14 days**, calculated as 24 hours of continuous operation per day.

The runtime projection uses the measured task times and the specified number of workers under continuous operation.

## 7. Practical implications

R = 5300 is set in run_config.R. The decision, its rule, the estimated planning MCSE, and the runtime are recorded in outputs/planning_summary/main_study_plan.csv. The other rounding options in the candidate table are comparison values.

Another complete precision pilot is not required for the decision made. The checks currently imply no necessary change to n_test = 5000, RF with 200 trees, or the SNR levels.

To evaluate the available results after extraction, execute in R:

``` r
source(file.choose())
```

Select start.R in the project folder. The script updates the summary and export without model fits. Alternatively, use source("update_planning.R") in the project directory.

All discussed pilot and check scripts are included. The next study phase is the actual main simulation with the agreed R and its confirmatory analysis.
