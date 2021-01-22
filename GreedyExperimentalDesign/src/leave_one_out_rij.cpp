#include <Rcpp.h>
using namespace Rcpp;
//using namespace arma;

/////// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector compute_avg_abs_rijss_left(NumericMatrix Rij) {
  int n_R = Rij.nrow();
  NumericVector avg_abs_rijss_left(n_R);
  
  int num_elements = (n_R - 1) * (n_R - 2);
  for (int i_left = 0; i_left < n_R; i_left++){
    double tot = 0;
    for (int i = 0; i < n_R; i++){
      for (int j = 0; j < n_R; j++){
        if (i != j && i != i_left && j != i_left){
          tot += Rij(i, j);
        }
      }
    }
    avg_abs_rijss_left[i_left] = tot / num_elements;
  }
  return avg_abs_rijss_left;
}
