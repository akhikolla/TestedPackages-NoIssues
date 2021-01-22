Title:       rpf

# Why use RPF?

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jpritikin/rpf.svg?branch=master)](https://travis-ci.org/jpritikin/rpf)
[![Codecov test coverage](https://codecov.io/gh/jpritikin/rpf/branch/master/graph/badge.svg)](https://codecov.io/gh/jpritikin/rpf?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/rpf?color=blue)](https://cran.r-project.org/package=rpf)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/rpf)](https://cranlogs.r-pkg.org/badges/rpf)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rpf)](https://cranlogs.r-pkg.org/badges/grand-total/rpf)
<!-- badges: end -->

The idea behind RPF is modularity. Most item factor analysis software
is not modular. Modularity facilitates more contributors and cross
pollination between projects.

# Installation

To get the current released version from CRAN:

```R
install.packages("rpf")
```

# Developer notes

There are a number of useful scripts in the `tools` subdir:

* `install` -- Installs the package as quickly as possible. Skips
  building the vignettes and documentation.

* `build` -- Builds a source tarball

* `check` -- Builds a source tarball and checks it

* `rox` -- Re-generates the documentation.

* `test` -- Runs the test suite using the uninstalled tests against the
  installed package.

* `autodep` -- Recalculates the header file dependences

If you're working on the C++ code, you probably want to adjust the
settings in `src/Makevars`.
