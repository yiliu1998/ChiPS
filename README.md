# Chips
`Chips` is an R package that implements the estimation of general causal estimands using the propensity score weighted (PSW) estimators. We employ two ways for the estimation:

- Model the propensity scores using the whole dataset, predict the propensity score values on all participants (estimate-then-plug-in), and calculate the PSW estimator. The variance estimation is by nonparametric bootstrap (see [Efron (1992)](https://link.springer.com/chapter/10.1007/978-1-4612-4380-9_41)).
- Employ the double/debiased machine learning (DML) algorithms by [Chernozhukov et al. (2018)](https://academic.oup.com/ectj/article/21/1/C1/5056401) for propensity score modelling on training sets and predicting on prediction set, using multiple sample-spltting and cross-fitting. This is called the DML-based PSW estimator in our paper. The variance estimation for this DML-based estimator is direct using the mean square of the influnce function. 

For models used for the propensity score, we allow conventional parametric models such as the logistic regression or other nonparametric and machine learning models through the `SuperLearner` package. Please see [Van der Laan et a. (2007)](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html) and [Polley et al. (2010)](https://biostats.bepress.com/ucbbiostat/paper266/?TB_iframe=true&width=370.8&height=658) for details about the `SuperLearner`. Using `SuperLearner::listWrappers()` in R can output all available options for modelling the propensity score. 

## Installation
To install the latest version of the R package from GitHub, please this code in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/Chips")
```

## Usage
The package contains two main functions `Chips_Boot()` and `Chips_DML()`. The input data of the function should be observational data with a binary treatment (`A`, valued from 0 and 1), an outcome variable (`Y`) and covariate matrix or data frame (`X`). Once the package is downloaded, please use the following code to view more detailed use of the package:

```r
library(Chips)
?Chips_Boot
?Chips_DML
```
As a remark, the current version of the package assumes **no missing data** in the input dataset. Therefore, if your data contain any missing values, we suggest either to conduct a complete-case analysis (e.g., if the missingness rate is small), or impute the missing data first using some existing methods such as multiple imputation. 

## Contact
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions).

## Reference
Please cite the following paper:
TBA
