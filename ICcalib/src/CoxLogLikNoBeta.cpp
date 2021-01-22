#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// beta- a value

//data is not assumed to be sorted by time

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double CoxLogLikNoBeta(arma::vec gamma, arma::vec tm, arma::vec event, arma::mat Z) {
  int n = tm.size();
  double denom=0;
  double logDenom=0;
  double logNumer=0;
  double logLik=0;
  int iCaseNum=-1;
  arma::vec GamZ = Z * gamma;
  arma::vec ExpGamZ = exp(Z * gamma);
  for (int i = 0; i < n; ++i)
  {
    if (event[i]) {
      iCaseNum += 1;
      logNumer += GamZ[i];
      denom = ExpGamZ[i];
      for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        denom += ExpGamZ[j];
        }
        }
     logDenom += log(denom);
    } }
  logLik = logNumer - logDenom; 
  
//NumericVector contrib=ps*exp(beta)+1-ps;
//List aa; aa["logLik"] = logLik;
//  return aa;
    return logLik;
  }


