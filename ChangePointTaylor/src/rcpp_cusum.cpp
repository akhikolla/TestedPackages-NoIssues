
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cusum(NumericVector x, double X_bar ) {
  int n = x.size();
  
  NumericVector S(n+1);
  S[0] =0;
  for(int i = 1; i <= n; ++i) {
    S[i] = S[i-1] + x[i-1]-X_bar;
  }
  return S;
}
