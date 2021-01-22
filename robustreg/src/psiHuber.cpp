#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector psiHuber_rcpp(NumericVector r, double c) {
int n =r.size();
NumericVector y = clone(r);

for(int i=0;i<n;i++){
 if(r[i]>c){
  y[i]=c;
 }
 else if(r[i]<-c){
   y[i]=-c;
 }
 else{
  y[i]=r[i];
 }
}
 return y;
}
 








