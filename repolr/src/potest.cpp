// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List potest(Rcpp::List mod, Rcpp::List hgmat, Rcpp::List X,
                 unsigned int categories, Rcpp::NumericVector ctimes){

    // inputs  
    arma::rowvec id = Rcpp::as<arma::rowvec>(mod["id"]);
    arma::rowvec y = Rcpp::as<arma::rowvec>(mod["y"]);
    arma::sp_mat Xmat = Rcpp::as<arma::sp_mat>(X["design"]);
    arma::rowvec beta = Rcpp::as<arma::rowvec>(mod["coefficients"]);
    arma::sp_mat irmat = Rcpp::as<arma::sp_mat>(hgmat["icormat"]);
    arma::rowvec times = Rcpp::as<arma::rowvec>(ctimes);
    unsigned int ntimes = times.n_elem;

    // calculate new variables
    arma::vec nexbeta = arma::vectorise(arma::repmat(beta.cols(categories - 1, beta.n_cols - 1), categories-1, 1));
    arma::mat exbeta = arma::join_rows(beta.cols(0, categories - 2), nexbeta.t());
    arma::mat linpred = Xmat * exbeta.t();
    arma::mat fitted = arma::exp(linpred) / (arma::exp(linpred) + 1);
    unsigned int maxid = arma::max(id);
    unsigned int nid = ntimes * (categories - 1);
    unsigned int nbeta = beta.n_cols;

    arma::mat varmat = arma::sqrt(fitted % (1 - fitted));
    arma::sp_mat dmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    arma::sp_mat vmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    vmat.diag() = 1 / varmat.t();
    dmat.diag() = arma::exp(linpred) / arma::pow(1 + arma::exp(linpred), 2.0);
    arma::sp_mat ddmat = Xmat.t() * dmat;
    
    // test statistic
    arma::mat ym = arma::conv_to<arma::mat>::from(y);
    arma::mat residuals = ym.t() - fitted;
    arma::sp_mat cprod1 = vmat * dmat * Xmat;
    arma::sp_mat cprod2 = irmat * vmat; 
    arma::mat umat = cprod1.t() * (cprod2 * residuals);
    arma::sp_mat vcovmat = vmat * irmat * vmat;
    arma::sp_mat wmat = ddmat * vcovmat * ddmat.t();
    arma::mat mwmat(wmat);
    double teststat = arma::as_scalar(umat.t() * arma::inv_sympd(mwmat)  * umat);
    unsigned int testdf = (categories - 2) * (nbeta - (categories - 1));



    // output
    return Rcpp::List::create(Rcpp::Named("umat")=umat,
                              Rcpp::Named("wmat")=wmat,
                              Rcpp::Named("teststat")=teststat,
                              Rcpp::Named("testdf")=testdf);

}


