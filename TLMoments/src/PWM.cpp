#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////
// COMPUTATION.METHOD = PWM

// C++ function to calculate probability-weighted moments
// [[Rcpp::export]]
double pwm_C(NumericVector x, int r) {
  int n = x.size();
  NumericVector xs(clone(x));
  xs = xs.sort();
  double sum = 0;
  double w;
  double vorfaktor = 1;

  for (int k = n; k >= (n-r); k--) {
    vorfaktor /= k;
  }
  for (int j = 1; j <= n; j++) {
    w = 1;
    for (int i = 1; i <= r; i++) {
      w *= (j-i);
    }
    sum += w * xs[j-1];
  }
  return sum * vorfaktor;
}
