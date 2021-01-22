
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/EMMIXgene)](https://cran.r-project.org/package=EMMIXgene)
[![Travis-CI Build
Status](https://travis-ci.org/andrewthomasjones/EMMIXgene_no_f.svg?branch=master)](https://travis-ci.org/andrewthomasjones/EMMIXgene_no_f)

# EMMIXgene

Provides unsupervised selection and clustering of microarray data using
mixture models. Following the methods described in McLachlan, Bean and
Peel (2002) <https://doi.org/10.1093/bioinformatics/18.3.413> a subset of genes are
selected based one the likelihood ratio statistic for the test of one
versus two components when fitting mixtures of t-distributions to the
expression data for each gene. The dimensionality of this gene subset is
further reduced through the use of mixtures of factor analyzers,
allowing the tissue samples to be clustered by fitting mixtures of
normal distributions.

## Installation

You can install EMMIXgene from github with:

``` r
# install.packages("devtools")
devtools::install_github("andrewthomasjones/EMMIXgene_no_f")
```

## Example
