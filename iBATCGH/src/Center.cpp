#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::mat Center(arma::mat Y){
int g=Y.n_cols;

	for (int i=0; i<g;i++){
   Y.col(i)=Y.col(i)-mean(Y.col(i));
}

return Y;

}
