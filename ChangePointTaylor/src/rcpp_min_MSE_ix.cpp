
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int min_MSE_ix(NumericVector x ) {
  int n = x.size();
  double pos_inf_proxy = std::numeric_limits<double>::max();
  
  NumericVector MSE_vals(n);
  MSE_vals[0] = pos_inf_proxy;
  MSE_vals[n-1] = pos_inf_proxy;
  NumericVector x_1;
  NumericVector x_2 ;
  
  for(int i =1; i < n-1; ++i) {
    x_1 = x[seq(0,i)];
    x_2 = x[seq(i+1,n-1)];
    MSE_vals[i] = sum(pow((x_1-mean(x_1)),2)) + sum(pow((x_2-mean(x_2)),2));
  }
  return which_min(MSE_vals)+1;
}
