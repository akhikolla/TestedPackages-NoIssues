// Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#include <Rcpp.h>
using namespace Rcpp;

double sumKernel(
    const NumericMatrix& X,    /* naip image */
  const NumericMatrix& W,    /* pre computed spatial weights */
  int i,      /* current location in rows */
  int j,      /* current location in columns */
  int nRow,   /* number of Rows */
  int nCol    /* number of Columns */
) {
  int dRow = W.nrow();
  int dCol = W.ncol();
  
  /* adjustment that must be applied for edge effects */
  int k, l, n;
  
  double sumYK = 0, sumK = 0;
  
  int k_local;
  int l_local;
  
  /* the starts */
  int k_start = std::max(i - dRow/2, (int)0);
  int l_start = std::max(j - dCol/2, (int)0);

  /* the stops */
  int k_stop = std::min(i + dRow/2 + 1, nRow);
  int l_stop = std::min(j + dCol/2 + 1, nCol);

  for(n = 0; n < X.nrow(); n++){
    for(k=k_start, k_local=k_start - i + (dRow/2); 
        k < k_stop; k++, k_local++) {
      for(l=l_start, l_local=l_start -j + (dCol/2);
          l < l_stop; l++, l_local++) {
        if(NumericVector::is_na(X(n, k * nCol + l))) continue;
        sumYK += X(n, k * nCol + l) * W(k_local, l_local);
        sumK += W(k_local, l_local);
      }
    }
  }
  if(sumK == 0.0){
    return NA_REAL;
  } else{
    return sumYK/sumK;
  }
}


// [[Rcpp::export]]
NumericVector mean_est(NumericMatrix X, int nRow, int nCol, NumericMatrix W) {
  // each row of X is a row stacked image
  // X.ncol() == nRow * nCol
  int i,j;
  NumericVector mu(nRow*nCol);
  
  for( i=0; i < nRow; i++) {
    for( j=0; j < nCol; j++) {
      Rcpp::Rcout << "i = " << i << std::endl;
      Rcpp::Rcout << "j = " << j << std::endl;
      mu[i*nCol + j] = sumKernel(X, W, i, j, nRow, nCol);
    }
  }
  return mu;
}
