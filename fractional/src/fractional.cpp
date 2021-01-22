#include <Rcpp.h>
using namespace Rcpp;

IntegerVector rat_one(double x, double eps, int maxConv) {
  int p0, p1 = 1, p2 = (int) floor(x),
      q0, q1 = 0, q2 = 1, b, i = 0;
  double z = x - (double) p2;
  while(++i < maxConv) {
    if(fabs(x - (double) p2 / (double) q2) < eps) break;
    z = 1/z; b = (int) floor(z); z = z - b;
    p0 = p1; p1 = p2; p2 = b*p1 + p0;
    q0 = q1; q1 = q2; q2 = b*q1 + q0;
  }
  return IntegerVector::create(p2, q2, i-1);
}

//' @describeIn ratr C++ version of the same function used for speed
//' @export
//' @import Rcpp
//' @useDynLib fractional
// [[Rcpp::export]]
IntegerMatrix rat(NumericVector x, double eps = 1.0e-6, int maxConv = 20) {
  int nx = x.length();
  IntegerMatrix PQC(nx, 3);
  PQC.attr("dimnames") = List::create(R_NilValue,
           CharacterVector::create("Pn", "Qn", "n"));
  for(int i = 0; i < nx; i++) {
    PQC(i, _) = rat_one(x[i], eps, maxConv);
  }
  return PQC;
}




