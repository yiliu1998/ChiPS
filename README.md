# GenCausalEst
R package implementing the estimation of general causal estimands in observational study using propensity score weighting estimators. The propensity score model in the current version of the package uses only the basic logsitic regression. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/General-Estimands")
```

## Use
The package contains only a main function `Gen_Estimand()` for use. The input data of the function should be observational data with a binary treatment (`A`, valued from 0 and 1), an outcome variable (`Y`) and covariate matrix or data frame (`X`). Once the package is downloaded, please use the following code to view more detailed use of the package:

```r
library(GenCausalEst)
?Gen_Estimand
```
In addition, please note that **we do not allow missing data** in the dataset for implementing our function. Therefore, if your data contain any missing values, we suggest either to conduct a complete-case analysis (e.g., if the missingness rate is small), or impute the missing data first. 

## Contact
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat[at]gmail[dot]com, if you have any questions).
