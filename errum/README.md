
<!-- README.md is generated from README.Rmd. Please edit that file -->

## errum

<!-- badges: start -->

[![R build
status](https://github.com/tmsalab/errum/workflows/R-CMD-check/badge.svg)](https://github.com/tmsalab/errum/actions)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN
status](https://www.r-pkg.org/badges/version/errum)](https://CRAN.R-project.org/package=errum)
<!-- badges: end -->

Perform a Bayesian estimation of the Exploratory reduced Reparameterized
Unified Model (‘ErRUM’) described by Culpepper and Chen (2018).

### Installation

You can install `errum` from CRAN using:

``` r
install.packages("errum")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("tmsalab/errum")
```

### Usage

To use the `errum` package, load it into *R* using:

``` r
library("errum")
```

From there, the errum model can be estimated using:

``` r
errum_model = errum(<data>, <k>, chain_length = 10000)
```

### Authors

James Joseph Balamuta, Steven Andrew Culpepper, and Jeffrey A. Douglas

### Citing the errum package

To ensure future development of the package, please cite `errum` package
if used during an analysis or simulation studies. Citation information
for the package may be acquired by using in *R*:

``` r
citation("errum")
```

### License

GPL (\>= 2)
