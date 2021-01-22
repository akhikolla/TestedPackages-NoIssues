#include <Rcpp.h>
using namespace Rcpp;

//' @keywords internal
//' Function to calculate distance between all training data points and one test point.
// [[Rcpp::export]]
NumericVector calc_distance_C(NumericMatrix trainX, NumericVector testX) {
  int C = trainX.ncol();
  int R = trainX.nrow();
  NumericVector dist(R);

  double temp = 0.;

  for(int i = 0; i < R; i++) {
    temp = 0.;
    for(int j = 0; j < C; j++) {
      temp += pow(trainX(i, j) - testX(j), 2.);
    }
    dist(i) = pow(temp, 1./2.);
  }

  return dist;
}
