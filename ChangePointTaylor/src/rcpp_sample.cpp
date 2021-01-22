
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sample_cpp(NumericVector x, int n) {
  return sample(x, n);
}
