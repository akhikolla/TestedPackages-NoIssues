// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List icormat(List mod, List smat, List cmat, Rcpp::String modtype){

    // inputs
    arma::sp_mat icmat = Rcpp::as<arma::sp_mat>(cmat["icmat"]);
    arma::sp_mat ismat = Rcpp::as<arma::sp_mat>(smat["ismat"]);
    arma::rowvec id;
    if (modtype == "glm"){
      Rcpp::List datalist = Rcpp::as<List>(mod["data"]);
      id = Rcpp::as<arma::rowvec>(datalist["subjects"]);
      } else if (modtype == "gee") {
      id = Rcpp::as<arma::rowvec>(mod["id"]);
     }  

    // calculate new variables
    unsigned int maxid = arma::max(id);
    unsigned int ntimes = icmat.n_rows;
    unsigned int cats = ismat.n_rows + 1;
    unsigned int nid = ntimes * (cats - 1);
    arma::mat mismat(ismat);
    arma::mat micmat(icmat);
     
    // block diagonal locations
    arma::mat clocblock = arma::linspace(0, nid * maxid-1, nid * maxid);
    arma::mat vrow = arma::repmat(clocblock,1,nid);
    clocblock.reshape(nid, maxid);
    arma::mat vcol = arma::repmat(clocblock.t(),1,nid);
    arma::vec loccol = arma::vectorise(vcol.t());
    arma::vec locrow = arma::vectorise(vrow.t());
    arma::mat locats = arma::join_rows(locrow,loccol);
    arma::umat ulocats = arma::conv_to<arma::umat>::from(locats.t());
   
    // calculate icormat
    arma::mat icor = arma::kron(micmat, mismat);
    arma::mat ficor = arma::repmat(icor,1, maxid);
    arma::vec vicor = arma::vectorise(ficor);
    arma::sp_mat icormat(ulocats, vicor);

     // output
    return Rcpp::List::create(Rcpp::Named("irmat")=icormat);

}

