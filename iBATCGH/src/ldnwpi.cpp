#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double ldnwpi(double x, double mean, double sigma) 
{

double out=-(pow((x-mean),2)/(2*pow(sigma,2)))-log(sigma);

return out;

}

