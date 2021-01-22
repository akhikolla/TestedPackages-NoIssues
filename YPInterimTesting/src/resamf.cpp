#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector awlrstat(int repnum, int n, NumericMatrix quan1, NumericMatrix quan2, NumericVector sd1, NumericVector sd2) {

  int nci = quan1.ncol();

  NumericMatrix out(repnum, nci);
  NumericVector w1_tilde(nci);
  NumericVector w2_tilde(nci);
  NumericVector aw1_tilde(nci);
  NumericVector aw2_tilde(nci);

  for (int i = 0; i < repnum; ++i) {

    NumericVector error = Rcpp::rnorm(n);

    for (int j = 0; j < nci; ++j) {
      for (int k = 0; k < n; ++k) {
        w1_tilde[j] += quan1(k, j)*error[k];
        w2_tilde[j] += quan2(k, j)*error[k];
      }
    }

    for (int j = 0; j < nci; ++j) {
        aw1_tilde[j] = std::abs(w1_tilde[j]/sd1[j]);
        aw2_tilde[j] = std::abs(w2_tilde[j]/sd2[j]);
        out(i, j) = std::max(aw1_tilde[j], aw2_tilde[j]);
        w1_tilde[j] = 0;
        w2_tilde[j] = 0;
    }

  }

  return(out);
}

