#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]
double pRedge(double gamma, double omega, int rgm, int rgma, double add){
double out=log(add*gamma+((rgm==rgma)?omega:0));
return out;
}

// [[Rcpp::interfaces(cpp)]]
double pRmiddle(double gamma, double omega1, double omega2, int rgm, int rgmm, int rgmp, double add){
double out=log(add*gamma+((rgm==rgmm)?omega1:0)+((rgm==rgmp)?omega2:0));
return out;
}

