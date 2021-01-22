#include "RcppArmadillo.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
//' Function to calculate dann distance between one training data points and one test point.
// [[Rcpp::export]]
arma::vec DANN_distance_C(const arma::rowvec & x0, const arma::rowvec & x1, const arma::mat &  sigma) {
  arma::vec distance = (x0- x1) * sigma * (x0 - x1).t();
  return(distance);
}
