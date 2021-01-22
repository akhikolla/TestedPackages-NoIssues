/*
 * check OpenMP availability and set desired number of threads
 */

#ifdef _OPENMP
#include <omp.h>
#endif

int openmp_threads = 1;

// [[Rcpp::export]]
DataFrame CPP_get_openmp_threads() {
  int num_threads = openmp_threads;
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 0;
#endif
  DataFrame res =
    DataFrame::create(_["available"] = max_threads > 0,
                      _["max"] = max_threads, 
                      _["threads"] = num_threads);
  res.attr("row.names") = "OpenMP";
  return res;
}

// [[Rcpp::export]]
void CPP_set_openmp_threads(int n) {
  if (n < 1) stop("internal error -- number of threads must be >= 1");
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n > max_threads) n = max_threads;
  openmp_threads = n;
#else
  if (n > 1) Rf_warning("OpenMP support not available");
  openmp_threads = 1;
#endif
}

