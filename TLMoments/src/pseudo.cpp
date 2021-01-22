#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pseudo_C(NumericVector x, NumericVector k) {
  int n = x.size();
  int K = k.size();

  NumericVector Emp(n);
  NumericMatrix vec1(n, K);
  NumericMatrix vec2(n, K);
  NumericMatrix out(n, K);

  for (int i = 0; i < n; i++) {
    Emp[i] = (double)(i+1)/n;
  }

  for (int j = 0; j < K; j++) {

    for (int i = 0; i < n; i++) {
      vec1(i, j) = x[i] * std::pow(Emp[i], k[j]);
    }

    vec2(n-1, j) = (double)1/n * x[n-1] * k[j] * std::pow(Emp[n-1], k[j]-1);
    for (int i = (n-2); i >= 0; i--) {
      vec2(i, j) = (double)1/n * x[i] * k[j] * std::pow(Emp[i], k[j]-1) + vec2(i+1, j);
    }

    for (int i = 0; i < n; i++) {
      out(i, j) = vec1(i, j) + vec2(i, j);
    }
  }

  return out;
}
