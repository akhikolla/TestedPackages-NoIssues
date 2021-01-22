#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>
#include <cmath>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::export]]
int pdist(const arma::mat &x, arma::mat &D){
  int n = x.n_rows;
  //int d = x.n_cols;
  D = arma::zeros(n,n);
  for(int i=0;i<n;++i){
    for(int j=i+1;j<n;++j){
      D(j,i) = D(i,j) = std::sqrt(arma::sum(arma::square(x.row(i)-x.row(j))));
      //D(j,i) = D(i,j) = std::sqrt(D(i,j));
    }
  }
  return(0);
}

//[[Rcpp::export]]
int Ucenter(arma::mat &D){
  int n = D.n_rows;
  arma::rowvec colmeanD = sum(D,0)/(n-2);
  double mD = sum(colmeanD)/(n-1);
  D.each_row() -= colmeanD;
  D.each_col() -= colmeanD.t();
  D+=mD;
  for(int i=0;i<n;++i) D(i,i) = 0;
  return(0);
}

//[[Rcpp::export]]
int Vcenter(arma::mat &D){
  arma::rowvec colmeanD = arma::mean(D,0);
  double mD = arma::mean(colmeanD);
  D.each_row() -= colmeanD;
  D.each_col() -= colmeanD.t();
  D+=mD;
  return(0);
}

//[[Rcpp::export]]
int centering(arma::mat& D,std::string type="V"){

  int n = D.n_rows,d1, d2;
  if(type=="U"){
    d1 = n-2; d2 = n-1;
  }else{
    d1 = n; d2 = n;
  }
  arma::rowvec colmeanD = sum(D,0)/d1;
  double mD = sum(colmeanD)/d2;
  D.each_row() -= colmeanD;
  D.each_col() -= colmeanD.t();
  D+=mD;
  if(type=="U") for(int i=0;i<n;++i) D(i,i) = 0;
  return(0);
}

//[[Rcpp::export]]
int centering_from_data(const arma::mat& x,arma::mat &D,std::string type="V"){

  pdist(x,D);
  centering(D,type);
  return 0;

}
