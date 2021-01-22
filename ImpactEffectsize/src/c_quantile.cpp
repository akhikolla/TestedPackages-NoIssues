#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_quantile(NumericVector x, NumericVector probs) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(probs - 0.000000001)];
}
