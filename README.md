# Chips
`Chips` is an R package that implements the estimation of general causal estimands in observational study using propensity score weighting estimators. We employ two ways for the causal effect estimations:

- Estimate the propensity scores using some posited models (we allow traditional models such as a logistic regression or nonparametric models through the SuperLearner) and plug it in to the 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/ChiPS")
```

## Usage
The package contains two main functions `ChiPS()` and `ChiPS_DML()`. The input data of the function should be observational data with a binary treatment (`A`, valued from 0 and 1), an outcome variable (`Y`) and covariate matrix or data frame (`X`). Once the package is downloaded, please use the following code to view more detailed use of the package:

```r
library(ChiPS)
?ChiPS
?ChiPS_DML
```
In addition, please note that **we do not allow missing data** in the dataset for implementing our function. Therefore, if your data contain any missing values, we suggest either to conduct a complete-case analysis (e.g., if the missingness rate is small), or impute the missing data first. 

## Contact
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions).
