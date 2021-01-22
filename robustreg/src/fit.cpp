// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat fit_rcpp(arma::mat X, arma::vec b) {
return(X*b);
}


