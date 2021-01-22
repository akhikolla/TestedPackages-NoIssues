// Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#include <Rcpp.h>
using namespace Rcpp;

int which_equal(NumericVector x, int y){
  // return idx of x: which x equals y
  int idx = -1;
  NumericVector::iterator it = std::find(x.begin(), x.end(), y);
  if (it != x.end())
    idx = std::distance(x.begin(), it);
  return idx;
}
bool contains(NumericVector x, int y) {
  int idx;
  idx = which_equal(x, y);
  return idx >= 0; 
}

std::vector<int> intersect(std::vector<int> x, NumericVector y){
  std::vector<int> v;
  // should make sure x and y are sorted; in our case it should be
  // sort(x.begin(), x.end());
  // sort(y.begin(), y.end());
  std::set_intersection(x.begin(),x.end(),y.begin(),y.end(), back_inserter(v));
  return v;
}

std::vector<int> nbr_(int ii, int nRow, int nCol, int nnr){
    // The "following" neighorhood pixel indexes (including ii)
    // (i, j) is the current point coordinates
    // bvec: to the Bottom of (i,j) vector
    // rvec: to the Right of (i,j) vector
    int k, l; //cursor coordinates
    std::vector<int> bvec, rvec, vec;
    int i = ii / nCol;
    int j = ii % nCol;
    // Rcpp::Rcout << "i = " << (int)i << std::endl;
    // Rcpp::Rcout << "j = " << (int)j << std::endl;
    for(l = j; l <= std::min(j + nnr, nCol-1); l++){
      rvec.push_back(i*nCol + l);
    }
    if(i < nRow - 1){
      for(k = i+1; k <= std::min(i + nnr, nRow-1); k++){
        for(l = std::max(j - nnr, (int)0); l <= std::min(j + nnr, nCol-1); l++){
          bvec.push_back(k*nCol + l);
        }
      }

      vec.reserve( rvec.size() + bvec.size()); // preallocate memory
      vec.insert( vec.end(), rvec.begin(), rvec.end() );
      vec.insert( vec.end(), bvec.begin(), bvec.end() );
      return vec;
    } else {
      return rvec;
    }
}

std::vector<int> nbr_(int ii, int nRow, int nCol, int dRow, int dCol){
  // All neiborhood pixel indexes (including ii) within a given dRow x dCol rectangle
  // where ii is the centroid of the rectangle.
  // nRow, nCol are the data matrix dimension
  // dRow, dCol are the Weight matrix dimension.
  // The output is sorted (small to large)
  int k, l;
  int i = ii / nCol;
  int j = ii % nCol;
  std::vector<int> vec;
  for(k = std::max(i - dRow/2, (int)0); k <= std::min(i + dRow/2, nRow-1); k++){
    for(l = std::max(j - dCol/2, (int)0); l <= std::min(j + dCol/2, nCol-1); l++){
      vec.push_back(k*nCol + l);
    }
  }
  return vec;
}

// [[Rcpp::export]]
IntegerVector nbr(int ii, int nRow, int nCol, int dRow, int dCol){
  return wrap(nbr_(ii, nRow, nCol, dRow, dCol));
}

double emp_cov_(const NumericMatrix &X, int ii, int jj){
  //sparse empirical covariance estimation
  double tmp, res = 0.0;
  int n = 0;
  for(int k = 0; k < X.nrow(); k++){
    tmp = X(k, ii) * X(k, jj);
    if(!NumericVector::is_na(tmp)){
      res += tmp;
      n += 1;
    }
  }
  if(n < 2)
    return NA_REAL;
  else
    return res/(n-1);
}

  
double lc_cov_(const NumericMatrix &X, const NumericMatrix &W,
               int ii, int jj, int nRow, int nCol){
  // sparse local constant covariance estimation for points i and j
  // X: data matrix, each row is a row stacked image
  // W: weight matrix
  // ii: first pixel index; jj: second pixel index
  int dRow = W.nrow();
  int dCol = W.ncol();

  double sumEEKK = 0.0, sumKK = 0.0;
  
  int k1, l1, k2, l2;
  int k1_local, l1_local, k2_local, l2_local;

  int i1 = ii / nCol; //row index for the first pixel
  int j1 = ii % nCol; //column index for the first pixel
  int i2 = jj / nCol;
  int j2 = jj % nCol;
  
  /* the starts */
  int k1_start = std::max(i1 - dRow/2, (int)0);
  int l1_start = std::max(j1 - dCol/2, (int)0);
  int k2_start = std::max(i2 - dRow/2, (int)0);
  int l2_start = std::max(j2 - dCol/2, (int)0);
  
  /* the stops */
  int k1_stop = std::min(i1 + dRow/2 + 1, nRow);
  int l1_stop = std::min(j1 + dCol/2 + 1, nCol);
  int k2_stop = std::min(i2 + dRow/2 + 1, nRow);
  int l2_stop = std::min(j2 + dCol/2 + 1, nCol);
  
  for(int n = 0; n < X.nrow(); n++){
    for(k1 = k1_start, k1_local = k1_start - i1 + (dRow/2); 
        k1 < k1_stop; k1++, k1_local++) {
      for(l1 = l1_start, l1_local=l1_start - j1 + (dCol/2);
          l1 < l1_stop; l1++, l1_local++) {
        if(NumericVector::is_na(X(n, k1 * nCol + l1)))
          continue;
        for(k2 = k2_start, k2_local = k2_start - i2 + (dRow/2); 
            k2 < k2_stop; k2++, k2_local++) {
          for(l2 = l2_start, l2_local=l2_start - j2 + (dCol/2);
              l2 < l2_stop; l2++, l2_local++) {
            if(NumericVector::is_na(X(n, k2 * nCol + l2)))
              continue;
            if(k1 == k2 && l1 == l2)
              continue;
            sumEEKK += X(n, k1 * nCol + l1) * W(k1_local, l1_local) * X(n, k2 * nCol + l2) * W(k2_local, l2_local);
            sumKK += W(k1_local, l1_local) * W(k2_local, l2_local);
          }
        }
      }
    }
  }
  if(sumKK == 0.0){ // a standalone point, in this case ii should be equal to jj
    for(int n = 0; n < X.nrow(); n++){
      sumEEKK = X(n, ii) * X(n, jj);
      sumKK += 1;
    }
    return sumEEKK/sumKK;
  } else{
    return sumEEKK/sumKK;
  }
}

// [[Rcpp::export]]
double lc_cov1_(const NumericMatrix &X, const NumericMatrix &W,
               int ii, int jj, int nRow, int nCol, NumericVector &pidx){
  // sparse local constant covariance estimation for points i and j
  // X: data matrix, each row is a row stacked image with black hole columns removed
  // W: weight matrix
  // ii: first pixel index; jj: second pixel index IN X
  int dRow = W.nrow();
  int dCol = W.ncol();
  
  double sumEEKK = 0.0, sumKK = 0.0;
  
  int k1, l1, k2, l2, idx1, idx2;
  int k1_local, l1_local, k2_local, l2_local;
  
  int i1 = (int)pidx[ii] / nCol; //row index for the first pixel
  int j1 = (int)pidx[ii] % nCol; //column index for the first pixel
  int i2 = (int)pidx[jj] / nCol;
  int j2 = (int)pidx[jj] % nCol;
  
  /* the starts */
  int k1_start = std::max(i1 - dRow/2, (int)0);
  int l1_start = std::max(j1 - dCol/2, (int)0);
  int k2_start = std::max(i2 - dRow/2, (int)0);
  int l2_start = std::max(j2 - dCol/2, (int)0);
  
  /* the stops */
  int k1_stop = std::min(i1 + dRow/2 + 1, nRow);
  int l1_stop = std::min(j1 + dCol/2 + 1, nCol);
  int k2_stop = std::min(i2 + dRow/2 + 1, nRow);
  int l2_stop = std::min(j2 + dCol/2 + 1, nCol);
  
  for(k1 = k1_start, k1_local = k1_start - i1 + (dRow/2); 
      k1 < k1_stop; k1++, k1_local++) {
    for(l1 = l1_start, l1_local=l1_start - j1 + (dCol/2);
        l1 < l1_stop; l1++, l1_local++) {
      idx1 = which_equal(pidx, k1 * nCol + l1);
      if(idx1 < 0)
        continue;
      for(k2 = k2_start, k2_local = k2_start - i2 + (dRow/2); 
          k2 < k2_stop; k2++, k2_local++) {
        for(l2 = l2_start, l2_local=l2_start - j2 + (dCol/2);
            l2 < l2_stop; l2++, l2_local++) {
          if(k1 == k2 && l1 == l2)
            continue;
          idx2 = which_equal(pidx, k2 * nCol + l2);
          if(idx2 < 0)
            continue;
          for(int n = 0; n < X.nrow(); n++){
            // if 'k1 * nCol + l1' is a point in black hole or it is a missing value, skip the current loop
            if(NumericVector::is_na(X(n, idx1)))
              continue;
            if(NumericVector::is_na(X(n, idx2)))
              continue;
            sumEEKK += X(n, idx1) * W(k1_local, l1_local) * X(n, idx2) * W(k2_local, l2_local);
            sumKK += W(k1_local, l1_local) * W(k2_local, l2_local);
          }
        }
      }
    }
  }
  
  if(sumKK == 0.0){
    // a standalone point, in this case ii should be equal to jj
    // in this case return the variance estimation
    for(int n = 0; n < X.nrow(); n++){
      if(NumericVector::is_na(X(n, ii)))
        continue;
      sumEEKK += X(n, ii) * X(n, jj);
      sumKK += 1;
    }
    return sumEEKK/(sumKK-1);
  } else{
    return sumEEKK/sumKK;
  }
}

// [[Rcpp::export]]
DataFrame sparse_emp_cov_est(NumericMatrix X, int nRow, int nCol, int nnr) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<int> ridx;
  std::vector<int> cidx;
  std::vector<int> ii_nbr;
  std::vector<double> value;

  int ii, jj;
  for(ii = 0; ii < nRow*nCol; ii++){
    ii_nbr = nbr_(ii, nRow, nCol, nnr);
    for(jj = 0; jj < ii_nbr.size(); jj++){
      ridx.push_back(ii + 1);
      cidx.push_back(ii_nbr[jj] + 1);
      value.push_back(emp_cov_(X, ii, ii_nbr[jj]));
    }
  }
  return DataFrame::create( 
    _["ridx"]  = ridx, 
    _["cidx"]  = cidx, 
    _["value"] = value
  );
}

// [[Rcpp::export]]
DataFrame sparse_emp_cov_est1(NumericMatrix X, int nRow, int nCol, int nnr, NumericVector pidx) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<int> ridx;
  std::vector<int> cidx;
  std::vector<int> ii_nbr;
  std::vector<double> value;

  int ii, jj, jj1;
  for(ii = 0; ii < pidx.size(); ii++){
    ii_nbr = nbr_(pidx[ii], nRow, nCol, nnr);
    ii_nbr = intersect(ii_nbr, pidx);
    for(jj = 0; jj < ii_nbr.size(); jj++){
      ridx.push_back(ii + 1);
      jj1 = which_equal(pidx, ii_nbr[jj]);
      cidx.push_back(jj1 + 1);
      value.push_back(emp_cov_(X, ii, jj1));
    }
  }
  return DataFrame::create( 
    _["ridx"]  = ridx, 
    _["cidx"]  = cidx, 
    _["value"] = value
  );
}

// [[Rcpp::export]]
DataFrame sparse_lc_cov_est(NumericMatrix X, NumericMatrix W, int nRow, int nCol, int nnr) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<int> ridx;
  std::vector<int> cidx;
  std::vector<int> ii_nbr;
  std::vector<double> value;
  
  int ii, jj;
  for(ii = 0; ii < nRow*nCol; ii++){
    ii_nbr = nbr_(ii, nRow, nCol, nnr); //pixel indexes following the ii-th pixel
    for(jj = 0; jj < ii_nbr.size(); jj++){
      ridx.push_back(ii + 1);
      cidx.push_back(ii_nbr[jj] + 1);
      value.push_back(lc_cov_(X, W, ii, ii_nbr[jj], nRow, nCol));
    }
  }
  return DataFrame::create( 
    _["ridx"]  = ridx, 
    _["cidx"]  = cidx, 
    _["value"] = value
  );
}

// [[Rcpp::export]]
DataFrame sparse_lc_cov_est1(NumericMatrix X, NumericMatrix W, int nRow, int nCol, int nnr, NumericVector pidx) {
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  std::vector<int> ridx;
  std::vector<int> cidx;
  std::vector<int> ii_nbr;
  std::vector<double> value;
  
  int ii, jj, jj1;
  for(ii = 0; ii < pidx.size(); ii++){
    ii_nbr = nbr_(pidx[ii], nRow, nCol, nnr);
    ii_nbr = intersect(ii_nbr, pidx);
      for(jj = 0; jj < ii_nbr.size(); jj++){
        ridx.push_back(ii + 1);
        jj1 = which_equal(pidx, ii_nbr[jj]);
        cidx.push_back(jj1 + 1);
        value.push_back(lc_cov1_(X, W, ii, jj1, nRow, nCol, pidx));
      }
  }
  return DataFrame::create( 
    _["ridx"]  = ridx, 
    _["cidx"]  = cidx, 
    _["value"] = value
  );
}
