#include <vector>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

/*
arma::mat pdist2(const arma::mat &x){
  int n = x.n_rows;
  arma::mat x2rs = repmat(arma::sum(arma::square(x),1),1,n);
  arma::mat D = x2rs+x2rs.t() - 2*x*x.t();
  D = clamp(D,+0.0,D.max()+1.0);
  return(D);
}
*/

arma::cube parccos(const arma::mat x){
  int n = x.n_rows;
  arma::cube a(n,n,n);
  for(int r=0;r<n;++r){
    arma::mat nx = x.each_row()-x.row(r);
    arma::colvec pdx = sqrt(sum(pow(nx,2),1));
    arma::mat pd_x = pdx*pdx.t();
    a.slice(r) = nx*nx.t();
    a.slice(r)/=pd_x;
  }
  a.replace(arma::datum::nan,1);
  a.replace(arma::datum::inf,1);
  a = arma::clamp(a,-1,1);
  a = arma::acos(a);
  return(a);
}

arma::cube Pcenter1(const arma::cube a){
  int n = a.n_rows;
  arma::mat a1(n,n), a2(n,n), b1(n,n), b2(n,n);
  arma::vec am(n), bm(n);
  a1 = mean(a,1);
  a2 = mean(a,0);
  for(int r=0;r<n;++r){
    am[r] = mean(mean(a.slice(r)));
  }
  arma::cube A=a;
  for(int i=0;i<n;++i){
    A.subcube(0,i,0,size(n,1,n))-=a1;
    A.subcube(i,0,0,size(1,n,n))-=a2;
    A.slice(i)+=am[i];
  }
  return(A);
}


arma::cube Pcenter(const arma::cube a){
  int n = a.n_rows;
  arma::mat a1(n,n), a2(n,n), b1(n,n), b2(n,n);
  arma::vec am(n), bm(n);
  a1 = mean(a,1);
  a2 = mean(a,0);
  for(int r=0;r<n;++r){
    am[r] = mean(mean(a.slice(r)));
  }
  arma::cube A=a;
  arma::mat A1(n,n),A2(n,n);
  for(int r=0;r<n;++r){
    A1.each_col() = a1.col(r);
    A2.each_row() = a2.col(r).t();
    A.slice(r)-=(A1+A2);
    A.slice(r)+=am[r];
  }
  return(A);
}


double parccov(const arma::cube A,const arma::cube B){
  int n = A.n_rows;
  arma::cube AB(n,n,n);
  AB = A%B;
  double pc=0;
  for(int r=0;r<n;++r){
    pc+=mean(mean(AB.slice(r)));
  }
  pc/=n;
  return(pc);
}

//[[Rcpp::export]]
double pcovCpp(const arma::mat x,const arma::mat y){
  //int n = x.n_rows;
  arma::cube a, b;
  a = parccos(x);
  b = parccos(y);
  arma::cube A, B;
  A = Pcenter(a);
  B = Pcenter(b);
  double pc = parccov(A,B);
  /*
  px = parccov(A,A);
  py = parccov(B,B);
  Rcpp::List result = Rcpp::List::create(Named("px")=px,
                                         Named("py")=py,
                                         Named("pxy")=pxy);
  */
  return(pc);
}


NumericVector mpcovCpp(const arma::mat &x,const arma::mat& y){

  int R = y.n_cols;
  arma::cube a,b;
  a = parccos(x);
  arma::cube A,B;
  A = Pcenter(a);
  NumericVector res(R);
  for(int i=0;i<R;++i){
    b = parccos(y.col(i));
    B = Pcenter(b);
    res[i] = parccov(A,B);
  }
  return(res);
}


List pcov_perm_Cpp(const arma::mat& x,const arma::mat& y,int R=500){
  int n = x.n_rows;
  arma::cube a, b;
  a = parccos(x);
  b = parccos(y);
  arma::cube A, B, AB(n,n,n),AA(n,n,n),BB(n,n,n);;
  A = Pcenter(a);
  B = Pcenter(b);
  double pc = parccov(A,B);
  arma::vec pc_perm(R);
  arma::uvec idx(n);
  for(int i=0;i<n;++i){idx[i]=i;}
  for(int i=0;i<R;++i){
    idx = arma::shuffle(idx);
    B = B.slices(idx);
    for(int j=0;j<n;++j){
      B.slice(j) = B.slice(j)(idx,idx);
    }
    pc_perm[i] = parccov(A,B);
  }
  double pval = 1.0-sum(pc>pc_perm)/(R+0.0);
  Rcpp::List result = Rcpp::List::create(Named("pc")=pc,
                                         Named("pc.perm")=pc_perm,
                                         Named("p.value")=pval);
  return(result);
}

