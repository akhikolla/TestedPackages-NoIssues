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


double sumdxy2(const arma::vec& x,const arma::vec& y){
  // x is not sorted and will not be sorted
  int n = x.size();
  arma::uvec idx = sort_index(x), idy = sort_index(y); //x[idx[i]] = i smallest value in x
  arma::uvec ranky(n);
  for(int i=0;i<n;++i) ranky[idy[i]] = i;
  double sumdxy = 0.0, sumy = arma::sum(y);
  arma::vec bit1 = arma::zeros(n+1), bit2 = arma::zeros(n+1);
  arma::vec rsm=arma::zeros(n), rsym=arma::zeros(n);
  for(int i=n-1;i>-1;--i){
    rsm[idx[i]] = psum(bit1,n-1-ranky[idx[i]]);
    rsym[idx[i]] = psum(bit2,n-1-ranky[idx[i]]);
    update(bit1,n,n-1-ranky[idx[i]],1);
    update(bit2,n,n-1-ranky[idx[i]],y[idx[i]]);
  }
  bit1.zeros(); bit2.zeros();
  double sm, sym;
  for(int i=0;i<n;++i){
    sm = 2*(psum(bit1,ranky[idx[i]]) + rsm[idx[i]]) - n;
    sym = 2*(psum(bit2,ranky[idx[i]]) + rsym[idx[i]]) - sumy;
    sumdxy += 2*x[idx[i]]*y[idx[i]]*sm;
    sumdxy -= 2*x[idx[i]]*sym;
    update(bit1,n,ranky[idx[i]],1);
    update(bit2,n,ranky[idx[i]],y[idx[i]]);

  }
  return(sumdxy);

}


arma::vec amean(const arma::vec& x){

  int n = x.size();
  arma::vec am(n);
  arma::uvec idx = sort_index(x);
  double sum = arma::sum(x), csum = 0;
  for(int i=0;i<n;++i){
    //rankx[idx[i]] = i;
    am[idx[i]] = x[idx[i]]*(i*2-n)/n - (2*csum-sum)/n;
    csum += x[idx[i]];
  }
  return am;

}

//[[Rcpp::export]]
double pdcov1v1v1(arma::vec x, arma::vec y, arma::vec z, std::string type){

  int n = x.n_elem;
  double d1,d2,d3;
  if(type=="U"){d1 = n*(n-3.0); d2 = (n-2)*(1-3.0/n); d3 = (1-1.0/n)*(1-2.0/n)*(1-3.0/n);
  }else{d1 = n*(n+0.0); d2 = n; d3 = 1;}
  double dxysum = sumdxy2(x,y);
  double dxzsum = sumdxy2(x,z);
  double dyzsum = sumdxy2(y,z);
  double sumz = sum(z);
  double dzzsum = 2*arma::sum(arma::square(z))*n - 2*sumz*sumz;
  arma::vec am = amean(x), bm = amean(y), cm = amean(z);
  double amm = arma::mean(am), bmm = arma::mean(bm), cmm = arma::mean(cm);
  double dxy = dxysum/d1 - 2*arma::sum(am%bm)/d2 + amm*bmm/d3;
  double dxz = dxzsum/d1 - 2*arma::sum(am%cm)/d2 + amm*cmm/d3;
  double dyz = dyzsum/d1 - 2*arma::sum(bm%cm)/d2 + bmm*cmm/d3;
  double dz2 = dzzsum/d1 - 2*arma::sum(cm%cm)/d2 + cmm*cmm/d3;
  double pdc = dxy - dxz*dyz/dz2;
  return pdc;

}

//[[Rcpp::export]]
double pdcor1v1v1(arma::vec x, arma::vec y, arma::vec z, std::string type){

  int n = x.n_elem;
  double d1,d2,d3;
  if(type=="U"){
    d1 = n*(n-3.0); d2 = (n-2)*(1-3.0/n); d3 = (1-1.0/n)*(1-2.0/n)*(1-3.0/n);
  }else{
    d1 = n*(n+0.0); d2 = n; d3 = 1;
  }
  double dxysum = sumdxy2(x,y);
  double dxzsum = sumdxy2(x,z);
  double dyzsum = sumdxy2(y,z);
  double sumx = arma::sum(x), sumy = arma::sum(y), sumz = arma::sum(z);
  double dxxsum = 2*arma::sum(arma::square(x))*n - 2*sumx*sumx;
  double dyysum = 2*arma::sum(arma::square(y))*n - 2*sumy*sumy;
  double dzzsum = 2*arma::sum(arma::square(z))*n - 2*sumz*sumz;
  arma::vec am = amean(x), bm = amean(y), cm = amean(z);
  double amm = arma::mean(am), bmm = arma::mean(bm), cmm = arma::mean(cm);
  double dxy = dxysum/d1 - 2*arma::sum(am%bm)/d2 + amm*bmm/d3;
  double dxz = dxzsum/d1 - 2*arma::sum(am%cm)/d2 + amm*cmm/d3;
  double dyz = dyzsum/d1 - 2*arma::sum(bm%cm)/d2 + bmm*cmm/d3;
  double dx2 = dxxsum/d1 - 2*arma::sum(am%am)/d2 + amm*amm/d3;
  double dy2 = dyysum/d1 - 2*arma::sum(bm%bm)/d2 + bmm*bmm/d3;
  double dz2 = dzzsum/d1 - 2*arma::sum(cm%cm)/d2 + cmm*cmm/d3;
  double pdxy = dxy - dxz*dyz/dz2;
  double pdxx = dx2 - dxz*dxz/dz2;
  double pdyy = dy2 - dyz*dyz/dz2;
  double pdr = pdxy/sqrt(pdxx*pdyy);
  if(isnan(pdr)||pdxx<0||pdyy<0) pdr = 0;
  return pdr;

}
