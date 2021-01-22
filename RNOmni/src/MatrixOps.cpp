// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Matrix Inner Product
//'
//' Calculates the product \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
// [[Rcpp::export]]
SEXP matIP(const arma::mat A, const arma::mat B){
  const arma::mat AtB = A.t()*B;
  return Rcpp::wrap(AtB);
}

//' Matrix Inverse
//' 
//' Calcualtes \eqn{A^{-1}}.
//'
//' @param A Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP matInv(const arma::mat A){
  const arma::mat Ai = arma::pinv(A);
  return Rcpp::wrap(Ai);
}

//' Schur complement
//'
//' Calculates the efficient information \eqn{I_{bb}-I_{ba}I_{aa}^{-1}I_{ab}}. 
//'
//' @param Ibb Information of target parameter
//' @param Iaa Information of nuisance parameter
//' @param Iba Cross information between target and nuisance parameters
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP SchurC(const arma::mat Ibb, const arma::mat Iaa,
            const arma::mat Iba){
  const arma::mat Ibba = Ibb-Iba*arma::solve(Iaa,Iba.t(),arma::solve_opts::likely_sympd);
  return Rcpp::wrap(Ibba);
}
