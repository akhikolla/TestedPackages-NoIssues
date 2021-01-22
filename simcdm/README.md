
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simcdm

[![Build
Status](https://travis-ci.org/tmsalab/simcdm.svg)](https://travis-ci.org/tmsalab/simcdm)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/simcdm)](https://cran.r-project.org/package=simcdm)
[![CRAN
Status](https://cranchecks.info/badges/worst/simcdm)](https://cran.r-project.org/web/checks/check_results_simcdm.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/simcdm?color=brightgreen)](http://www.r-pkg.org/pkg/simcdm)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/simcdm?color=brightgreen)](http://www.r-pkg.org/pkg/simcdm)
[![Coverage
status](https://codecov.io/gh/tmsalab/simcdm/branch/master/graph/badge.svg)](https://codecov.io/github/tmsalab/simcdm?branch=master)

The goal of `simcdm` is to provide flexible ways to simulate data under
cognitive diagnostic models.

## Installation

You can install simcdm from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("tmsalab/simcdm")
```

## Usage

To use `simcdm`, load the package using:

``` r
library("simcdm")
```

## Overview

There are four distinct sets of functions within the package:

  - Attributes: `attribute_classes()`, `attribute_bijection()`,
    `attribute_inv_bijection()`, and `sim_subject_attributes()`.
  - Matrix: `sim_q_matrix()` and `sim_eta_matrix()`
  - Deterministic Input, Noisy And Gate (DINA): `sim_dina_items()` and
    `sim_dina_attributes()`
  - reduced Reparameterized Unified Model (rRUM): `sim_rrum_items()`

Functions that use random numbers to simulate values are named with the
prefix of `sim_*()`. This is done to allow for functions to be quickly
identified and used through autocomplete inside of the [RStudio
IDE](https://www.rstudio.com). At a later time, the `attribute_*()` will
likely be moved to a different package.

For more details, please see the [package vignettes](vignettes/):

  - [Overview of `simcdm`](vignettes/overview-simcdm.html)
  - [Using `simcdm` in R packages](vignettes/simcdm-in-packages.html)

## Authors

James Joseph Balamuta and Steven Andrew Culpepper with contributions
from Aaron Hudson.

## Citing the `simcdm` package

To ensure future development of the package, please cite `simcdm`
package if used during the analysis or simulations. Citation information
for the package may be acquired by using in *R*:

``` r
citation("simcdm")
```

## License

GPL (\>= 2)
