#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]    
NumericVector psiBS_rcpp(NumericVector r, double c) {
  int n = r.size();
  NumericVector y = clone(r);
  
  for (int i=0; i<n; i++) {
    if (abs(r[i]) <=c ) {
      y[i] = (r[i]*pow(1-pow(r[i]/c,2),2));
    }
    else{
      y[i]=0;
    }
  }
  return y;
}




