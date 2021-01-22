# HRM 1.2.0

[![CRANstatus](https://www.r-pkg.org/badges/version/HRM)](https://cran.r-project.org/package=HRM)
[![](https://cranlogs.r-pkg.org/badges/HRM)](https://cran.r-project.org/package=HRM)
[![Travis-CI Build Status](https://travis-ci.org/happma/HRM.svg?branch=test)](https://travis-ci.org/happma/HRM)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/happma/HRM?branch=test&svg=true)](https://ci.appveyor.com/project/happma/HRM)
[![codecov](https://codecov.io/gh/happma/HRM/branch/test/graph/badge.svg)](https://codecov.io/gh/happma/HRM)


R package for analysing high-dimensional repeated measures for factorial designs. A description of this package can be found in [1], theoretical derivations of the test statistics are in [2] and [3].



To install the current development version:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("happma/HRM", ref = "test", dependencies = TRUE)
library(HRM)
```

With this package it is possible to test for main and interaction effects of up to three whole- or subplot-factors. In total, a maximum of four factors can be used. There are two different S3 methods available. The first method requires a list of matrices in the wide table format. The second method requires a data.frame in the long table format.

``` r
## hrm_test with a list of matrices

# number patients per group
n = c(10,10)
# number of groups
a=2
# number of variables
d=40

# defining the list consisting of the samples from each group
mu_1 = mu_2 = rep(0,d)
# autoregressive covariance matrix
sigma_1 = diag(d)
for(k in 1:d) for(l in 1:d) sigma_1[k,l] = 1/(1-0.5^2)*0.5^(abs(k-l))
sigma_2 = 1.5*sigma_1
X = list(mvrnorm(n[1],mu_1, sigma_1), mvrnorm(n[2],mu_2, sigma_2))
X=lapply(X, as.matrix)

hrm_test(data=X, alpha=0.05)


## hrm.test with a data.frame using a 'formula' object

# using the EEG dataset
hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)
```

To get confidence intervals for each factor combination you can use the generic function 'confint' for an object of class 'HRM'. This function calculates simultaneous confidence intervals which maintains the family wise error rate (FWER).
See the following code:

``` r
# using the EEG dataset
z <- hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)

# calculate 99% confidence intervals
confint(z, level = 0.99)

```

In the data there are 4 variables with each 10 regions. We can use a multivariate approach as the variables are on different scales. For that, we can use the function 'hrm_test' with the argument 'variable' set to the column name which contains the factor variable for the variables.

``` r
# using the EEG dataset
hrm_test(value ~ group*region, subject = subject, variable = variable, data = EEG)
```


Additionally, the package can be used with a GUI.
``` r
hrm_GUI()
```

## References

[1] Happ, M., Harrar, S. W., and Bathke, A. C. (2018). HRM: An R Package for Analysing High-dimensional Multi-factor Repeated Measures. The R Journal 10(1), 534--548. <a href="https://journal.r-project.org/archive/2018/RJ-2018-032/index.html">https://journal.r-project.org/archive/2018/RJ-2018-032/index.html</a>

[2] Happ, M., Harrar S. W. and Bathke, A. C. (2017). High-dimensional Repeated
  Measures. Journal of Statistical Theory and Practice. 11(3), 468-477. URL:
  <a href="https://doi.org/10.1080/15598608.2017.1307792">doi:10.1080/15598608.2017.1307792</a>.
  
[3] Happ, M., Harrar, S. W., & Bathke, A. C. (2016). Inference for low‐and high‐dimensional multigroup repeated measures designs with unequal covariance matrices. Biometrical Journal, 58(4), 810-830. <a href = "https://doi.org/10.1002/bimj.201500064">doi:10.1002/bimj.201500064</a>
