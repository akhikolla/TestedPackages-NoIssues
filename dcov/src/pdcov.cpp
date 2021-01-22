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
double pdcov(const arma::mat &x,const arma::mat &y,const arma::mat &z,std::string type="V"){

  if(x.n_cols==1&&y.n_cols==1&&z.n_cols==1){
    return pdcov1v1v1(x,y,z,type);
  }

  int n = x.n_rows;
  //int p = x.n_cols, q = y.n_cols;
  arma::mat Dx(n,n),Dy(n,n),Dz(n,n);
  centering_from_data(x,Dx,type);
  centering_from_data(y,Dy,type);
  centering_from_data(z,Dz,type);
  Dx = Dx - arma::sum(arma::sum(Dx%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  Dy = Dy - arma::sum(arma::sum(Dy%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  int d = type=="V"?(n*(n+0.0)):(n*(n-3.0));
  return arma::sum(arma::sum(Dx%Dy))/d;
}

//[[Rcpp::export]]
double pdcor(const arma::mat &x,const arma::mat &y,const arma::mat &z,std::string type="V"){

  if(x.n_cols==1&&y.n_cols==1&&z.n_cols==1){
    return pdcor1v1v1(x,y,z,type);
  }

  int n = x.n_rows;
  //int p = x.n_cols, q = y.n_cols;
  arma::mat Dx(n,n),Dy(n,n),Dz(n,n);
  centering_from_data(x,Dx,type);
  centering_from_data(y,Dy,type);
  centering_from_data(z,Dz,type);
  Dx = Dx - arma::sum(arma::sum(Dx%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  Dy = Dy - arma::sum(arma::sum(Dy%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  int d = type=="V"?(n*(n+0.0)):(n*(n-3.0));
  double pdc = sum(sum(Dx%Dy))/d;
  double pdx = sum(sum(Dx%Dx))/d;
  double pdy = sum(sum(Dy%Dy))/d;
  double pdr = pdc/std::sqrt(pdx*pdy);
  if(isnan(pdr)||pdx<0||pdy<0){
    return 0;
  }else return pdr;
}

