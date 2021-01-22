//' @useDynLib RTCC, .registration = TRUE
//' @importFrom Rcpp sourceCpp

#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dist1 (const NumericVector x){
  int n = x.size();
  float total = 0;
  for (int i = 0; i < n ; ++i) {
    for (int j = i + 1; j < n; ++j){
      total += fabs(x(i)-x(j));
    }
  }
  return total;
}
