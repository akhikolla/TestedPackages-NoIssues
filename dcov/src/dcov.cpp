#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include "dcov.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::export]]
double dcov(const arma::mat &x,const arma::mat &y, std::string type="V"){

  int n = x.n_rows, p = x.n_cols, q = y.n_cols;
  if(p==1&&q==1){
    return dcov1v1(x,y,type);
  }else{
    int d = type=="V"?(n*(n+0.0)):(n*(n-3.0));
    arma::mat Dx(n,n), Dy(n,n);
    centering_from_data(x,Dx,type);
    centering_from_data(y,Dy,type);
    double dcov = arma::sum(arma::sum(Dx%Dy))/d;
    return dcov;
  }
}

//[[Rcpp::export]]
double dcor(const arma::mat &x,const arma::mat &y, std::string type){

  int n = x.n_rows, p = x.n_cols, q = y.n_cols;
  if(p==1&&q==1){
    return dcor1v1(x,y,type);
  }else{
    int d = type=="V"?(n*(n+0.0)):(n*(n-3.0));
    arma::mat Dx(n,n), Dy(n,n);
    centering_from_data(x,Dx,type);
    centering_from_data(y,Dy,type);
    double dcov = arma::sum(arma::sum(Dx%Dy))/d;
    double dx2  = arma::sum(arma::sum(Dx%Dx))/d;
    double dy2  = arma::sum(arma::sum(Dy%Dy))/d;
    double dcor = dcov/sqrt(dx2*dy2);
    if(isnan(dcor)||dx2<0||dy2<0){
      return 0;
    }else return dcor;
  }

}

