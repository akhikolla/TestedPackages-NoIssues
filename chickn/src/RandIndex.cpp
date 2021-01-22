// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#include <RcppParallel.h>
#include <RcppArmadillo.h>

#include <iostream>

using namespace arma;
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//' @title RandIndex
//' @description Rand Index computation
//' @param Data is a two columns matrix. First column contains the clustering assignment result, 
//' second column ground truth clustering assignment each row correspond to the data vector index.
//' @return The rand index value
//' @export
//' @keywords internal
// [[Rcpp::export]]
List RandIndex(NumericMatrix Data){
  
size_t n = Data.nrow();
 double a=0, b=0, d = 0, c = 0;
  
  for(size_t i = 0; i< n; i++) {
    for (size_t j = 0; j < i; j++) {
      if((Data(i,0)==Data(j,0)) & (Data(i,1)==Data(j,1))){
         a++;
      }else if((Data(i,0)!=Data(j,0)) & (Data(i,1)!=Data(j, 1))){
          b++;
        
      }else if((Data(i,0)==Data(j,0)) & (Data(i,1)!=Data(j, 1))){
        c++;
      }else{
        d++; 
      }
    }
  }
  return List::create(
    _["TP"] = a,
    _["TN"] = b,
    _["FP"] = c,
    _["FN"] = d,
    _["RI"] = (a+b)/(a+b+c+d)
  );
}
