#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
NumericVector loopC(NumericVector nodelist,int al,IntegerVector x, IntegerVector x2, List pm_lc, NumericMatrix Lix, int finind){
  NumericVector::iterator k;
  LogicalVector res;
  IntegerVector daughter;
  int temp1, temp2;
  int k2;
  int l = 0;
  int g = 0;
  int n = x.size();
  LogicalVector in1(x2.size());
  LogicalVector in2(x2.size());
  IntegerVector indices = seq_len(n); //from 1 to n...
  for(k = nodelist.begin(); k != nodelist.end(); ++k) {
    k2  = nodelist[l];
    res = x == k2;
    daughter = x2[res];
    in1 = x2==daughter[0];
    in2 = x2==daughter[1];
    temp1 = as<int>(indices[in1]);
    temp2 = as<int>(indices[in2]);
    NumericMatrix pmtmp1 = pm_lc[temp1 - 1];
    NumericMatrix pmtmp2 = pm_lc[temp2 - 1];
    for(g = 0; g<al; ++g){
      Lix(k2 - 1, g) = sum(Lix.row(daughter[0] - 1) * pmtmp1.row(g)) *
      sum(Lix.row(daughter[1] - 1) * pmtmp2.row(g));
    }
    l=l+1;
  }
  NumericVector res2 = Lix(finind,_);
  return res2;
}
// [[Rcpp::export]]
NumericVector colMeansC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(ncol);
  
  for (int i = 0; i < ncol; i++) {
    double total = 0;
    for (int j = 0; j < nrow; j++) {
      total += x(j, i);
    }
    out[i] = total/nrow;
  }
  return out;
}
// [[Rcpp::export]]
NumericVector cat_loopC(int nocat, NumericVector nodelist,int al,IntegerVector x, IntegerVector x2, List pm_clc, NumericMatrix Lix, int finind){
  NumericMatrix res3(nocat,al);
  NumericVector ret(al);
  //  NumericVector retu(al);
  int h = 0;
  for(h = 0; h<nocat; ++h){
    res3(h,_) = loopC(nodelist, al, x, x2, pm_clc[h], Lix, finind);
  }
  ret = colMeansC(res3);
  //  retu = log(sum(ret*rootpr));
  return ret;
}
// [[Rcpp::export]]
NumericMatrix part_loopC(int nocat, NumericVector nodelist,int al,IntegerVector x, IntegerVector x2, List pm_plc, List Lix, int finind){
  int k = 0;
  int lenpat = Lix.size();
  NumericMatrix mat_plc(lenpat,al); 
  for(k = 0; k<lenpat; ++k){
    NumericMatrix Lixi = Lix[k];
    mat_plc(k,_) = cat_loopC(nocat, nodelist, al, x, x2, pm_plc, Lixi, finind);
  }
  return mat_plc;
}
// [[Rcpp::export]]
arma::mat expm(const arma::mat& M) {
  arma::mat  N = arma::expmat(M);
  return N;
}
