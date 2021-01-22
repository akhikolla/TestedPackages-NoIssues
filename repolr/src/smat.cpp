// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

List smat(Rcpp::NumericVector coeff){

    // input data
    arma::rowvec mcoeff(coeff);
    int categ = mcoeff.size();
    arma::sp_mat ismat(categ, categ);

    // construct matrices
    arma::rowvec pexpcoeff = arma::exp(mcoeff);
    arma::rowvec nexpcoeff = arma::exp(-mcoeff);
    arma::mat rs_mat = pexpcoeff.t() * nexpcoeff;
    arma::mat s_mat = arma::symmatu(arma::sqrt(rs_mat));
    arma::mat is_mat = arma::inv_sympd(s_mat);

    // construct sparse matrices
    arma::vec vs_mat = arma::vectorise(s_mat, 0);
    arma::umat rseq = arma::repmat(arma::linspace<arma::umat>(0, categ-1, categ), categ, 1);
    arma::umat cseq = arma::sort(rseq);
    arma::umat local = arma::join_rows(rseq,cseq);
    arma::sp_mat smat(local.t(), vs_mat);
    ismat.diag(0) = is_mat.diag(0);
    ismat.diag(1) = is_mat.diag(1);
    ismat.diag(-1) = is_mat.diag(-1);

    // output data
    return Rcpp::List::create(Rcpp::Named("smat")=smat,
                              Rcpp::Named("ismat")=ismat);

}


