<img src="inst/extdata/deploid.png" width="180">

[![License (GPL version 3)](https://img.shields.io/badge/license-GPL%20version%203-brightgreen.svg)](http://opensource.org/licenses/GPL-3.0)
[![Build Status](https://travis-ci.org/DEploid-dev/DEploid-r.svg?branch=master)](https://travis-ci.org/DEploid-dev/DEploid-r)
[![Build Status](https://ci.appveyor.com/api/projects/status/hi1nq97d5l68qs4r?svg=true)](https://ci.appveyor.com/project/shajoezhu/deploid-r)
[![Coverage Status](https://coveralls.io/repos/github/DEploid-dev/DEploid-r/badge.svg?branch=master)](https://coveralls.io/github/DEploid-dev/DEploid-r?branch=master)
[![codecov](https://codecov.io/gh/DEploid-dev/DEploid-r/branch/master/graph/badge.svg)](https://codecov.io/gh/DEploid-dev/DEploid-r)
[![CRAN RStudio Mirror Downloads](http://cranlogs.r-pkg.org/badges/DEploid)](https://cran.r-project.org/package=DEploid)

DEploid R package -- Deconvolute Mixed Genomes with Unknown Proportions
=================

Traditional ‘phasing’ programs are limited to diploid organisms. Our method modifies Li and Stephen’s algorithm with Markov chain Monte Carlo (MCMC) approaches, and builds a generic framework that allows haloptype searches in a multiple infection setting.


Installation
------------

It is recommended to use the current CRAN version. It can be installed from within R using

```R
> install.packages('DEploid')
```

For the developing version, please install Rcpp package first. From the R-console, type

```R
> install.packages("Rcpp")
```

(NOTE: If you are using Windows, please install `Rtools` from [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/))

then

```R
> install.packages("devtools")
> library(devtools)
> install_github("mcveanlab/DEploid-r")
```

Usage
-----

Please refer to the help page and examples for each function. For example,
```R
> library(DEploid)
> ?dEploid
> ?plotProportions
```

Licence
-------

You can freely use all code in this project under the conditions of the GNU GPL Version 3 or later.


Citation
--------

Citation
--------

If you use `dEploid` with the flag `-ibd`, please cite the following paper:

Zhu, J. S., J. A. Hendry, J. Almagro-Garcia, R. D. Pearson, R. Amato, A. Miles, D. J. Weiss, T. C. D. Lucas, M. Nguyen, P. W. Gething, D. Kwiatkowski, G. McVean, and for the Pf3k Project. (2018) The origins and relatedness structure of mixed infections vary with local prevalence of *P. falciparum* malaria. *eLife*, 40845, doi: https://doi.org/10.7554/eLife.40845.


If you use `dEploid` in your work, please cite the program:

Zhu, J. S., J. A. Garcia, G. McVean. (2018) Deconvolution of multiple infections in *Plasmodium falciparum* from high throughput sequencing data. *Bioinformatics* 34(1), 9-15. doi: https://doi.org/10.1093/bioinformatics/btx530.
