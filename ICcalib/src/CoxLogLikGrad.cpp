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

arma::vec CoxLogLikGrad(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Z) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
  
  arma::vec Szero =  arma::zeros(sumD);
  arma::mat Sone = arma::zeros(sumD,nPars);
  int iCaseNum=-1;
  int jCaseNum=-1;
  double beta = theta[0];
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::mat contrib=1 + ps*(exp(beta)-1);
  arma::vec GamZ = Z * gamma;
  arma::vec ExpGamZ = exp(Z * gamma);
  
  arma::mat nu = 1 + ps*(exp(beta)-1);
  arma::mat nuDerivBeta= ps*exp(beta);
  
  
  arma::vec grad = arma::zeros(nPars);
  
  // First a loop to calculate s0 and s1
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
     iCaseNum += 1;
     Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamZ[i];
     Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,i)*ExpGamZ[i];
     for(int iGam = 1; iGam < nPars; ++iGam) {
       Sone(iCaseNum,iGam) += nu(iCaseNum,i)*ExpGamZ[i]*Z(i,iGam-1);
      }
     for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamZ[j];
        Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,j)*ExpGamZ[j];
        for(int jGam = 1; jGam < nPars; ++jGam) {
          Sone(iCaseNum,jGam) += nu(iCaseNum,j)*ExpGamZ[j]*Z(j,jGam-1);
        }
       }
      }
    }
  }
  
   for (int k = 0; k < n; ++k) {
     if (event[k]) {
       jCaseNum += 1;
       grad[0] += nuDerivBeta(jCaseNum,k)/nu(jCaseNum,k) - Sone(jCaseNum,0)/Szero[jCaseNum];
       for (int kGam = 1; kGam < nPars; ++kGam) {
         grad[kGam] += Z(k,kGam-1) - Sone(jCaseNum,kGam)/Szero[jCaseNum];
       }
     }
   }
  
    return grad;
  }


