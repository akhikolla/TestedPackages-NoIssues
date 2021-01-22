#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
arma::rowvec dir(arma::rowvec a)
{

int l=a.n_cols;
arma::rowvec x(l);

GetRNGstate();
for (int i=0;i<l;i++){
x[i]=Rf_rgamma(a[i],1);
}
PutRNGstate();

double sm=sum(x);

return x/sm;
}

