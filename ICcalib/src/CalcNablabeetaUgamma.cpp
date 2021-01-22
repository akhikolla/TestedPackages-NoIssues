#include <RcppArmadillo.h>
using namespace Rcpp;

// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// psderiv - A cube (R array maybe) the last dimention is the eta paramter so the first slice is the matrix of derivatives with 
// respect to eta_1, for all people and all event times.
// beta- a value
//data is not assumed to be sorted by time
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CalcNablabeetaUgamma(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Z, arma::mat psDeriv) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
//  int nEta = psDeriv.n_slices;
  int iCaseNum=-1;
  int jCaseNum=-1;
  // double FirstTerm=0;
  // double SecondTerm=0;
  // double FirstSumType=0;
  // double SecondSumType=0;
  // double ThirdSumType=0;
  double beta = theta[0];
  
  //arma::vec all=0;
  //double all=0;
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
  arma::vec a =  arma::zeros(sumD);
  arma::vec ExpGamZ = exp(Z * gamma);
  arma::vec ExpGamZbeta = exp(Z * gamma + beta);
  arma::mat deriv = arma::zeros(1,nPars-1);
  arma::mat SecondTermNumer = arma::zeros(1,nPars-1);
  arma::mat SzeroEta = arma::zeros(sumD,nPars-1);
  arma::mat SoneEta = arma::zeros(sumD,nPars-1);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);


  // First a loop to calculate all sums (e.g., S0 and S1) that do not involve psDeriv
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
      iCaseNum += 1;
      Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamZ[i];
      arma::mat Zi = Z(i,arma::span::all);
      SzeroEta(iCaseNum, arma::span::all) +=   nu(iCaseNum,i)*ExpGamZbeta[i];
     SoneEta(iCaseNum, arma::span::all) +=   psDeriv(iCaseNum,i)*ExpGamZ[i]*Zi;
  StwoGammaBeta(iCaseNum,arma::span::all) +=  Zi*ExpGamZbeta[i]*psDeriv(iCaseNum,i);
  for(int j = 0; j < n; ++j) {
    if (tm[j]>tm[i]) {
      Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamZ[j];
      arma::mat Zj = Z(j,arma::span::all);
    SzeroEta(iCaseNum, arma::span::all) +=   nu(iCaseNum,j)*ExpGamZbeta[j];
    SoneEta(iCaseNum, arma::span::all) +=   psDeriv(iCaseNum,j)*ExpGamZ[j]*Zj;
      StwoGammaBeta(iCaseNum,arma::span::all) +=  Zj*ExpGamZbeta[j]*psDeriv(iCaseNum,j);
    }
   }
    }
  }

 for (int k = 0; k < n; ++k)
   {
     if (event[k]) {
       jCaseNum += 1;
       SecondTermNumer(0, arma::span::all) = SoneEta(jCaseNum, arma::span::all)%SzeroEta(jCaseNum, arma::span::all);
       deriv(0, arma::span::all) -= (StwoGammaBeta(jCaseNum,arma::span::all)*Szero[jCaseNum]- SecondTermNumer)/(Szero[jCaseNum]*Szero[jCaseNum]);
       }
     }
     return deriv;
   }

