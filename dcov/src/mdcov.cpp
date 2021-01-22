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
arma::rowvec mdcov(const arma::mat y, const arma::mat x, std::string type="V"){

  int n = x.n_rows, p = x.n_cols, q = y.n_cols;
  arma::rowvec dc(p);
  if(q>1){
    int d = type=="V"?(n*n):(n*(n-3));
    arma::mat Dx(n,n), Dy(n,n);
    centering_from_data(y,Dy,type);
    for(int j=0;j<p;++j){
      centering_from_data(x.col(j),Dx,type);
      dc[j] = arma::sum(arma::sum(Dx%Dy))/d;
    }
  }else{
    arma::vec yv = y;
    for(int j=0;j<p;++j){
      dc[j] = dcov1v1(x.col(j),yv,type);
    }
  }
  return dc;
}


//[[Rcpp::export]]
arma::rowvec mdcor(const arma::mat y, const arma::mat x, std::string type="V"){

  int n = x.n_rows, p = x.n_cols, q = y.n_cols;
  arma::rowvec dr(p);
  if(q>1){
    int d = type=="V"?(n*n):(n*(n-3));
    double dcov,dx2,dy2;
    arma::mat Dx(n,n), Dy(n,n);
    centering_from_data(y,Dy,type);
    dy2 = arma::sum(arma::sum(Dy%Dy))/d;
    for(int j=0;j<p;++j){
      centering_from_data(x.col(j),Dx,type);
      dcov = arma::sum(arma::sum(Dx%Dy))/d;
      dx2 = arma::sum(arma::sum(Dx%Dx))/d;
      dr[j] = dcov/std::sqrt(dx2*dy2);
      if(isnan(dr[j])||dx2<0||dy2<0){
        dr[j]=0;
      }
    }
  }else{
    arma::vec yv = y;
    for(int j=0;j<p;++j){
      dr[j] = dcor1v1(x.col(j),yv,type);
    }
  }
  return dr;
}

