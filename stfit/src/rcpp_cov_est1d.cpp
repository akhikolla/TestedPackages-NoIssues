// Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#include <Rcpp.h>
using namespace Rcpp;

double vecmin(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  return *it;
}

double vecmax(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  return *it;
}

//' Local constant covariance estimation
//' @param ids a vector indicating subject/group ids
//' @param time integer vector of observed time points, the minimum time unit is 1
//' @param resid vector of residual values used for covariance calculation
//' @param W weight vector, it contains both kernel and bandwidth information in general 
//' local polynomial estimation setting up
//' @param t1 time point 1
//' @param t2 time point 2
//' @retrun covariance value between t1 and t2
//' @export
// [[Rcpp::export]]
double lc_cov_1d(const NumericVector &ids, const NumericVector &time, const NumericVector &resid, 
                  const NumericVector &W, int t1, int t2){
  int W_size = W.size();
  double sumEEKK = 0.0, sumKK = 0.0;
  int N = ids.size();
  int time_min = (int)vecmin(time);
  int time_max = (int)vecmax(time);

  /* the starts */
  int k1_start = std::max(t1 - W_size/2, time_min);
  int k2_start = std::max(t2 - W_size/2, time_min);
  
  /* the stops */
  int k1_stop = std::min(t1 + W_size/2 + 1, time_max);
  int k2_stop = std::min(t2 + W_size/2 + 1, time_max);

  for(int i = 0; i < N; i++){
    if((time[i] >= k1_start) & (time[i] < k1_stop)){
      for(int j = 0; j < N; j++){
        if(i == j)
          continue;
        if(ids[i] == ids[j]){
          if((time[j] >= k2_start) & (time[j] < k2_stop)){
            sumEEKK += resid[i]*resid[j]*W[time[i] - t1 + W_size/2]*W[time[j] - t2 + W_size/2];
            sumKK += W[time[i] - t1 + W_size/2]*W[time[j] - t2 + W_size/2];
          }
        }
      }
    }
  }
  if(sumKK == 0.0){ // no points within the bandwidth
    Rcpp::Rcout << "sumKK is 0" << std::endl;
    return NA_REAL;
  } else{
    return sumEEKK/sumKK;
  }
}

//' Local constant covariance estimation
//' @param ids a vector indicating subject/group ids
//' @param time integer vector of observed time points, the minimum time unit is 1
//' @param resid vector of residual values used for covariance calculation
//' @param W weight vector, it contains both kernel and bandwidth information in general 
//' local polynomial estimation setting up
//' @param tt time vector
//' @retrun a covariance matrix evaluated at time points \code{tt} on the covariance function 
//' @export
// [[Rcpp::export]]
NumericMatrix lc_cov_1d_est(const NumericVector &ids, const NumericVector &time, const NumericVector &resid, 
                  const NumericVector &W, const NumericVector &tt){
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  int tt_size = tt.size();
  NumericMatrix out(tt_size, tt_size);

  for(int i = 0; i < tt_size; i++){
    // if(i == 0)
    //   Rcpp::Rcout << "tt[i]: " << tt[i] << std::endl;
    for(int j = 0; j <= i; j++){
       out(i,j) = lc_cov_1d(ids, time, resid, W, tt[i], tt[j]);
      if(j < i)
        out(j,i) = out(i,j);
      // if(NumericVector::is_na(out(i,j)))
      //   Rcpp::Rcout << "i: " << i << "; j: " << j<< std::endl;
    }
  }
  return(out);
}
