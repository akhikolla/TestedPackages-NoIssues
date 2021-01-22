#include <RcppArmadillo.h>
using namespace Rcpp;


//
// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// theta- a vector of coefficients, the first is beta

//data is not assumed to be sorted by time

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat CoxLogLikHess(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Z) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
  int iCaseNum=-1;
  int jCaseNum=-1;
  
  double beta = theta[0];
  double NablaBetaUbeta=0;
  double DenomNablaGammaUgamma = 0;
  
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
  arma::vec a =  arma::zeros(sumD);
  arma::vec GamZ = Z * gamma;
  arma::vec ExpGamZ = exp(Z * gamma);
  
  arma::mat Sone = arma::zeros(sumD,nPars);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat contrib=1 + ps*(exp(beta)-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);
  arma::mat nuDerivBeta= ps*exp(beta);
  arma::mat nua = 1 - ps;
  arma::mat FirstBetaTermNumer = (nuDerivBeta % (1-ps));
  arma::mat FirstBetaTermDenom = nu%nu;
  arma::mat FirstBetaTerm = FirstBetaTermNumer/FirstBetaTermDenom;
  arma::mat NablaGammaUgamma = arma::zeros(nPars-1,nPars-1);
  arma::mat NumerNablaGammaUgamma = arma::zeros(nPars-1,nPars-1);
  arma::mat NumerNablaGammaUbeta  = arma::zeros(nPars-1,1);
  arma::mat NablaGammaUbeta  = arma::zeros(1,nPars-1);
  arma::mat Hess = arma::zeros(nPars,nPars);
  
  arma::cube StwoGamma = arma::zeros(nPars-1,nPars-1,sumD);
  // First a loop to calculate all sums (e.g., S1 and S1)
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
     iCaseNum += 1;
     Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamZ[i];
     Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,i)*ExpGamZ[i];
     a[iCaseNum] += nua(iCaseNum,i)*ExpGamZ[i];
     arma::mat Zi = Z(i,arma::span::all);
     StwoGammaBeta(iCaseNum,arma::span::all) +=  Zi*ExpGamZ[i]*nuDerivBeta(iCaseNum,i);
     StwoGamma.slice(iCaseNum) += Zi.t() * Zi * Szero[iCaseNum];
     Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,i)*ExpGamZ[i]*Zi;
     for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamZ[j];
        Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,j)*ExpGamZ[j];
        a[iCaseNum] += nua(iCaseNum,j)*ExpGamZ[j];
        arma::mat Zj = Z(j,arma::span::all);
        StwoGammaBeta(iCaseNum,arma::span::all) +=  Zj*ExpGamZ[j]*nuDerivBeta(iCaseNum,j);
        StwoGamma.slice(iCaseNum) += Zj.t() * Zj * nu(iCaseNum,j)*ExpGamZ[j];
        Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,j)*ExpGamZ[j]*Zj;
       }
      }
    }
  }
  for (int k = 0; k < n; ++k) {
     if (event[k]) {
       jCaseNum += 1;
       NablaBetaUbeta += FirstBetaTerm(jCaseNum,k) - (Sone(jCaseNum,0)*a[jCaseNum])/(Szero[jCaseNum]*Szero[jCaseNum]);
       arma::mat SoneTerm = Sone(jCaseNum,arma::span(1,nPars-1));
       NumerNablaGammaUgamma = Szero[jCaseNum]*StwoGamma.slice(jCaseNum) -   SoneTerm.t()*SoneTerm;
       DenomNablaGammaUgamma = Szero[jCaseNum]*Szero[jCaseNum];
       NablaGammaUgamma -= NumerNablaGammaUgamma/DenomNablaGammaUgamma;
       NumerNablaGammaUbeta = Szero[jCaseNum]*StwoGammaBeta(jCaseNum,arma::span::all) - SoneTerm*Sone(jCaseNum,0);
       NablaGammaUbeta -= NumerNablaGammaUbeta/DenomNablaGammaUgamma;
       }         
  }
  Hess(0,0) = NablaBetaUbeta;
  Hess(arma::span(1,nPars-1),arma::span(1,nPars-1)) = NablaGammaUgamma; 
  Hess(0,arma::span(1,nPars-1)) = NablaGammaUbeta;
  Hess(arma::span(1,nPars-1),0) = NablaGammaUbeta.t();
  return Hess;
}


