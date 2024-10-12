# ChiPS
`ChiPS` is an R package that implements the estimation of general causal estimands using the propensity score weighted (PSW) estimator. 

To implement the PSW estimator, we first model the propensity score using the whole dataset, predict its values on all participants (estimate-then-plug-in). The variance estimation is by nonparametric bootstrap (see [Efron (1992)](https://link.springer.com/chapter/10.1007/978-1-4612-4380-9_41)). For the model used for the propensity score, we allow conventional parametric models such as the logistic regression or other nonparametric and machine learning models through the `SuperLearner` package. Please see [Van der Laan et a. (2007)](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html) and [Polley et al. (2010)](https://biostats.bepress.com/ucbbiostat/paper266/?TB_iframe=true&width=370.8&height=658) for details about the `SuperLearner`. Using `SuperLearner::listWrappers()` in R can output all available options for modelling the propensity score. 

The use of nonparametric machine learning models and methods ensemble in `SuperLearner` can be time-consuming and computationally heavy in your experiments. As such, we recommend users to employ parallel computation in their simulation studies. Please visit [this page](https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf) for helpful information about how to conduct parallel computing in R for saving running time. 

Our package also output the plots of propensity score distributions by treatment group, covariate balancing (the loveplot by absolute standardized mean difference [ASMD]), and effective sample size (ESS), which are standard practices in causal inference for checking covariate balance and examine the method efficiency. 

## Installation
To install the latest version of the R package from GitHub, please this code in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/ChiPS")
```

## Usage
The package contains a main function `ChiPS()`. The input data of the function should be observational data with a binary treatment (`A`, valued from 0 and 1), an outcome variable (`Y`) and covariate matrix or data frame (`X`). Once the package is downloaded, please use the following code to view more detailed use of the package:

```r
library(ChiPS)
?ChiPS
```
As a remark, the current version of the package assumes **no missing data** in the input dataset. Therefore, if your data contain any missing values, we suggest either to conduct a complete-case analysis (e.g., if the missingness rate is small), or impute the missing data first using some existing methods such as multiple imputation. 

## Contact
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions).

## Reference
Please cite the following paper:
TBA
