#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double truncnorm(double mean,double sd,double lower,double upper){

GetRNGstate();
double res= Rf_qnorm5(Rf_pnorm5(lower,mean,sd,1,0)+unif_rand()*(Rf_pnorm5(upper,mean,sd,1,0)-Rf_pnorm5(lower,mean,sd,1,0)),0,1,1,0);
PutRNGstate();

return ((res*sd)+mean);       
       
}

