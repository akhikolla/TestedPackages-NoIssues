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
Rcpp::List dcor_test(const arma::mat& x, const arma::mat& y, int R, std::string type="V"){

  int n = x.n_rows, p = x.n_cols, q = y.n_cols;
  double stat;
  arma::rowvec replicates(R);
  if(p==1&&q==1){
    stat = dcor1v1(x,y,type);
    for(int i=0;i<R;++i){
      replicates[i] = dcor(x,shuffle(y),type);
    }
  }else{
    arma::mat Dx(n,n);
    arma::mat Dy(n,n);
    centering_from_data(x,Dx,type);
    centering_from_data(y,Dy,type);
    double dx2 = arma::sum(arma::sum(Dx%Dx));
    double dy2 = arma::sum(arma::sum(Dy%Dy));
    double dxy = std::sqrt(dx2*dy2);
    stat = arma::sum(arma::sum(Dx%Dy))/dxy;
    if(isnan(stat)||dx2<0||dy2<0){ stat = 0; }
    arma::uvec idx(n);
    for(int i=0;i<n;++i){idx[i]=i;}
    for(int i=0;i<R;++i){
      idx = arma::shuffle(idx);
      Dy = Dy(idx,idx);
      replicates[i] = arma::sum(arma::sum(Dx%Dy))/dxy;
      if(isnan(replicates[i])){
        replicates[i] = i%2?1e-16:-1e-16;
      }
    }
  }
  double p_val = (1.0+arma::sum(stat<replicates))/(1.0+R);
  Rcpp::List res = Rcpp::List::create(Named("statistic")=stat,
                                      Named("replicates")=replicates,
                                      Named("p.values")=p_val);
  return res;
}

//[[Rcpp::export]]
Rcpp::List pdcor_test(const arma::mat& x, const arma::mat& y, const arma::mat& z, int R, std::string type="U"){

  int n = x.n_rows;
  //int p = x.n_cols, q = y.n_cols;
  arma::mat Dx(n,n),Dy(n,n),Dz(n,n);
  centering_from_data(x,Dx,type);
  centering_from_data(y,Dy,type);
  centering_from_data(z,Dz,type);
  Dx = Dx - arma::sum(arma::sum(Dx%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  Dy = Dy - arma::sum(arma::sum(Dy%Dz))/arma::sum(arma::sum(Dz%Dz))*Dz;
  int d = type=="V"?n*n:n*(n-3);
  double pdc = sum(sum(Dx%Dy))/d;
  double pdx = sum(sum(Dx%Dx))/d;
  double pdy = sum(sum(Dy%Dy))/d;
  double pdxy = std::sqrt(pdx*pdy);
  double stat = pdc/pdxy;
  if(isnan(stat)||pdx<0||pdy<0){
    stat = 0;
  }
  arma::rowvec replicates(R);
  arma::uvec idx(n);
  for(int i=0;i<n;++i){idx[i]=i;}
  for(int i=0;i<R;++i){
    idx = arma::shuffle(idx);
    Dy = Dy(idx,idx);
    replicates[i] = arma::sum(arma::sum(Dx%Dy))/d/pdxy;
    if(isnan(replicates[i])){
      replicates[i] = i%2?1e-16:-1e-16;
    }
  }
  double p_val = (1.0+arma::sum(stat<replicates))/(1.0+R);
  Rcpp::List res = Rcpp::List::create(Named("statistic")=stat,
                                      Named("replicates")=replicates,
                                      Named("p.values")=p_val);
  return res;
}
