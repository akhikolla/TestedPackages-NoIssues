#ifndef DCUTILS_H
#define DCUTILS_H

# include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat expVec(double x, int deg);

arma::mat cMat (int k, NumericVector phi);

arma::mat invBMat (int k);

#endif
