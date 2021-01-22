
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double S_diff(NumericVector S) {
  return max(S) - min(S);
}
