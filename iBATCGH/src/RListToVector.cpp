#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::vec RListToVector(Rcpp::List xList,int G, int T){

int n = xList.size();
arma::vec vector(G*T);
vector.zeros();

for(int i=0; i<n; i++){
	vector(Rcpp::as<arma::uvec>(xList[i]))=vector(Rcpp::as<arma::uvec>(xList[i]))+1;
}

return vector;
}
