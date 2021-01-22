#include <RcppArmadillo.h>

// [[Rcpp::export]]
double computeMARCH(arma::mat e, int h){
  
  int K = e.n_cols, N = e.n_rows, pos = 0, eeVechCols = (K + 1 ) * K / 2;
  double tempValue = 0;
  
  arma::mat eeVech(N, eeVechCols, arma::fill::zeros);
  arma::mat eeVechDemeaned, resi;
  arma::mat regressors(N, 1 + h * eeVechCols, arma::fill::zeros);
  
  regressors(arma::span::all, 0).fill(1);
  
  
  for(int i = 0; i < N - h; ++i){
    pos = 0;
    for(int j = 0; j < K; ++j){
      for(int k = j; k < K; ++k){
        
        tempValue = e(i, j) * e(i, k);
        eeVech(i, pos) = tempValue;
        
        for(int l = 0; l < h; ++l){
          regressors(i + l + 1, l * eeVechCols + pos + 1) = tempValue;
        }
        pos++;
      }
    }
  }
  
  for(int i = N - h; i < N; ++i){
    pos = 0;
    for(int j = 0; j < K; ++j){
      for(int k = j; k < K; ++k){
        
        tempValue = e(i, j) * e(i, k);
        eeVech(i, pos) = tempValue;
        
        for(int l = 0; l < h; ++l){
          if(i + l + 1 < N) regressors(i + l + 1, l * eeVechCols + pos + 1) = tempValue;
        }
        pos++;
      }
    }
  }
  
  resi = eeVech - regressors * arma::inv(regressors.t() * regressors) * regressors.t() * eeVech;
  
  eeVechDemeaned = eeVech;
  
  for(size_t i = 0; i < eeVechDemeaned.n_cols; ++i) eeVechDemeaned.col(i) = eeVechDemeaned.col(i) - mean(eeVechDemeaned.col(i));
  
  
  return .5 * N * K * (K + 1) - N * arma::trace(resi.t() * resi * inv(eeVechDemeaned.t() * eeVechDemeaned));
}
