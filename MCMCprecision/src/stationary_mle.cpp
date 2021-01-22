#include <Rcpp.h>
using namespace Rcpp;


// estimation under assumption that MC is reversible
//
// Trendelkamp-Schroer, B., Wu, H., Paul, F., & Noé, F. (2015).
// Estimation and uncertainty of reversible Markov models.
// The Journal of Chemical Physics, 143(17), 174101.
// https://doi.org/10.1063/1.4934536
// [[Rcpp::export]]
NumericVector stationary_reversible (NumericVector pi, NumericMatrix N,
                                     double abstol = 1e-5, int maxit = 1e5)
{
  NumericVector pi0(clone(pi));
  int M = N.cols();
  NumericVector N_row = rowSums(N);

  int cnt = 0;
  double diff = 1.;
  while ((diff > abstol) && (cnt < maxit))
  {
    pi0 = clone(pi);
    for (int i = 0; i < M; i++)
    {
      pi[i] = sum( (N.row(i) + N.column(i))/(N_row[i]/pi0[i] + N_row/pi0));
    }
    cnt += 1;
    diff = max(abs(pi0 - pi));
    // Rcout << cnt << " / ", diff << "::: " << pi << " // " << pi0 << "\n";
  }
  if (cnt == maxit)
    warning("Maximum number of iterations reached.");
  return (pi);
}
