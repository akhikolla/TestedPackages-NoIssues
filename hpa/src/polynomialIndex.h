#ifndef hpa_polynomialIndex_H
#define hpa_polynomialIndex_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

NumericMatrix polynomialIndex(NumericVector pol_degrees,
                              bool is_validation);

String printPolynomial(NumericVector pol_degrees, 
								       NumericVector pol_coefficients,
								       bool is_validation);

#endif
