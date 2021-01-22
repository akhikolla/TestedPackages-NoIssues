# Spbsampling

An [R](https://www.r-project.org) package for *spatially balanced
sampling*.

The **Spbsampling** package provides functions to draw *spatially
balanced samples*. In particular, the implemented sampling designs allow
to select probability samples well spread over the population of
interest, in any dimension and using any distance function
(e.g. Euclidean distance, Manhattan distance).

For details regarding the implemented sampling designs, you may want to
look at the references section.

The implementation has been done in C++ through the use of **Rcpp** and
**RcppArmadillo**.

Authors: Francesco Pantalone, Roberto Benedetti, Federica Piersimoni.

Maintainer: Francesco Pantalone.

## Installation

You can install the released version of Spbsampling from
[CRAN](https://CRAN.R-project.org)

``` r
install.packages("Spbsampling")
```

## References

Benedetti R, Piersimoni F (2017). A spatially balanced design with
probability function proportional to the within sample distance.
*Biometrical Journal*, 59(5), 1067–1084.
<https://doi.org/10.1002/bimj.201600194>

Benedetti R, Piersimoni F (2017). Fast Selection of Spatially Balanced
Samples. *arXiv*. <https://arxiv.org/abs/1710.09116>

## License

GPL-3
