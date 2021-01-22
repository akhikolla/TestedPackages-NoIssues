// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List ordgee(Rcpp::List mod, Rcpp::List icormat, Rcpp::List X, 
                  Rcpp::String modtype, Rcpp::NumericVector ctimes, unsigned int categories,
                  unsigned int omaxit, double otol){

    // inputs
    arma::sp_mat Xmat = Rcpp::as<arma::sp_mat>(X["design"]);
    arma::rowvec linpred = Rcpp::as<arma::rowvec>(mod["linear.predictors"]);
    arma::rowvec fitted = Rcpp::as<arma::rowvec>(mod["fitted.values"]);
    arma::rowvec y = Rcpp::as<arma::rowvec>(mod["y"]);
    arma::rowvec residuals = y - fitted;
    arma::rowvec beta = Rcpp::as<arma::rowvec>(mod["coefficients"]);
    arma::rowvec times = Rcpp::as<arma::rowvec>(ctimes);
    unsigned int ntimes = times.n_elem;
    arma::rowvec id;
    if (modtype == "glm"){
      Rcpp::List datalist = Rcpp::as<List>(mod["data"]);
      id = Rcpp::as<arma::rowvec>(datalist["subjects"]);
      } else if (modtype == "gee") {
      id = Rcpp::as<arma::rowvec>(mod["id"]);
     }   
      
    // calculate new variables
    unsigned int maxid = arma::max(id);
    unsigned int ncoeff = beta.n_elem;
    unsigned int nid = ntimes * (categories - 1);
    arma::sp_mat irmat(nid * maxid, nid * maxid);
    if (modtype == "glm"){
      irmat = arma::speye<arma::sp_mat>(nid * maxid, nid * maxid);
      } else if (modtype == "gee") {
      irmat = Rcpp::as<arma::sp_mat>(icormat["irmat"]);
    }  
    arma::rowvec varmat = arma::sqrt(fitted % (1 - fitted));
    arma::sp_mat vmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    arma::sp_mat dmat = arma::zeros<arma::sp_mat>(nid * maxid, nid * maxid);
    vmat.diag() = 1 / varmat;
    dmat.diag() = arma::exp(linpred) / arma::pow(1 + arma::exp(linpred), 2.0);
    arma::mat nbeta = arma::conv_to<arma::mat>::from(beta);
    arma::mat ilinpred = arma::conv_to<arma::mat>::from(linpred);
    arma::mat iy = arma::conv_to<arma::mat>::from(y);
    arma::mat ifitted = arma::conv_to<arma::mat>::from(fitted);
    arma::mat iresiduals = arma::conv_to<arma::mat>::from(residuals);
    arma::sp_mat cprod1(nid * maxid, ncoeff), cprod2(nid * maxid, nid * maxid);
    arma::sp_mat ivcovmat(ncoeff, ncoeff);
    arma::mat eqgee(ncoeff, 1), updatebeta(ncoeff, 1);


    // loop until stop loop is true
    arma::mat crit, bloop;
    bool convergence = false, stoploop = false;
    unsigned int count = 0;
    arma::mat ibeta = nbeta;

      while(stoploop == false){
      
       // update beta
       count = count + 1;   
       for(unsigned int i=1; i<=omaxit; i++){
         ilinpred = Xmat * nbeta.t();
         dmat.diag() = arma::exp(ilinpred) / arma::pow(1 + arma::exp(ilinpred), 2.0);
         ifitted = arma::exp(ilinpred) / (arma::exp(ilinpred) + 1);
         iresiduals = iy - ifitted.t();
         vmat.diag() = 1 / arma::sqrt((ifitted % (1 - ifitted)));
         cprod1 = vmat * dmat * Xmat;
         cprod2 = irmat * vmat; 
         ivcovmat = cprod1.t() * (cprod2 * dmat * Xmat);    
         eqgee = cprod1.t() * (cprod2 * iresiduals.t());
         // superlu may be more efficient for large beta
         updatebeta = arma::spsolve(ivcovmat, eqgee, "lapack");
         bloop = arma::abs(updatebeta.t() / nbeta);
         if( bloop.max() < 100*otol){break;}
         nbeta = nbeta + updatebeta.t();

       }

      // control loop and while
      crit = arma::abs((ibeta-nbeta)/nbeta);
      if(crit.max() < otol){convergence = true; stoploop = true;}
      if(count >= omaxit){stoploop = true;}
      ibeta = nbeta;

      }

      // recalculate data for output
      ilinpred = Xmat * nbeta.t();
      ifitted = arma::exp(ilinpred) / (arma::exp(ilinpred) + 1);
      linpred = arma::vectorise(ilinpred,1);
      fitted = arma::vectorise(ifitted,1);

      // output
      return Rcpp::List::create(Rcpp::Named("y")=y,
                                Rcpp::Named("fitted.values")=fitted,
                                Rcpp::Named("linear.predictors")=linpred,
                                Rcpp::Named("id")=id,
                                Rcpp::Named("max.id")=nid,
                                Rcpp::Named("coefficients")=nbeta,
                                Rcpp::Named("convergence")=convergence);

}


