#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::colvec initial(arma::mat A) 
{

int n=A.n_rows;
arma::rowvec x(n);
x.ones();
arma::mat eye(n,n);
eye.eye();

arma::mat ones(n,n);
ones.ones();

arma::colvec output(n);
output=arma::trans((x)*arma::inv(eye-A+ones));

return output;
}

