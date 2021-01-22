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
double CalcNablabeetaUbeta(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Z, arma::mat psDeriv) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
//  int nEta = psDeriv.n_slices;
  int iCaseNum=-1;
  int kCaseNum=-1;
  double FirstTerm=0;
  double SecondTerm=0;
//  double FirstSumType=0;
//  double SecondSumType=0;
//  double ThirdSumType=0;
  double beta = theta[0];
  
  //arma::vec all=0;
  double all=0;
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
  arma::vec a =  arma::zeros(sumD);
  arma::vec ExpGamZ = exp(Z * gamma);
  arma::vec ExpGamZbeta = exp(Z * gamma + beta);
  
  arma::mat Sone = arma::zeros(sumD,nPars);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);
  arma::mat nua = 1 - ps;
  arma::mat nuaExpbeta = nua*exp(beta);
  
  arma::mat nuDerivBeta= ps*exp(beta);
  
  // First a loop to calculate all sums (e.g., S1 and S1) that do not involve psDeriv
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
      iCaseNum += 1;
      Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamZ[i];
      a[iCaseNum] += nua(iCaseNum,i)*ExpGamZ[i];
      for(int j = 0; j < n; ++j) {
        if (tm[j]>tm[i]) {
          Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamZ[j];
          a[iCaseNum] += nua(iCaseNum,j)*ExpGamZ[j];
        }
      }
    }
  }
  
  
//  for (int iEta = 0; iEta < nEta-1; ++iEta) {
    arma::mat SderivEta = arma::zeros(sumD);
   // arma::mat psDerivEta = psDeriv.slice(iEta);
    arma::mat psDerivEta = psDeriv;
    arma::mat TermOneNumerDeriv= nuaExpbeta%psDerivEta;
  // A a loop to calculate all the sums that do involve psDeriv
  int jCaseNum=-1;
//  int kCaseNum=-1;
    for (int j = 0; j < n; ++j) {
    if (event[j]) {
      jCaseNum += 1;
      SderivEta[jCaseNum] +=  psDerivEta(jCaseNum,j)*ExpGamZbeta[j];
      for(int k = 0; k < n; ++k) {
        if (tm[k]>tm[j]) {
          SderivEta[jCaseNum] +=  psDerivEta(jCaseNum,k)*ExpGamZbeta[k];
        }
      }
    }
  }


  for (int l = 0; l < n; ++l)
  {
    if (event[l]) {
      kCaseNum += 1;
      FirstTerm += TermOneNumerDeriv(kCaseNum,l)/(nu(kCaseNum,l)*nu(kCaseNum,l));
      SecondTerm += (SderivEta[kCaseNum]-a[kCaseNum])/(Szero[kCaseNum]*Szero[kCaseNum]);
      //(1-ps(iCaseNum,i))*contribDerivPbeta(iCaseNum,i)/(nu(iCaseNum,i)*nu(iCaseNum,i));
      // FirstSumType = nuDerivBeta(iCaseNum,i);
      // SecondSumType = nu(iCaseNum,i);
      // ThirdSumType = contribDerivPbeta(iCaseNum,i);
      // for(int j = 0; j < n; ++j) {
      //  if (tm[j]>tm[i]) {
      //    FirstSumType += nuDerivBeta(iCaseNum,j);
      //    SecondSumType += nu(iCaseNum,j);
      //    ThirdSumType += contribDerivPbeta(iCaseNum,j);
      //   }
      //   }
      //  SecondTerm += (ThirdSumType*SecondSumType-FirstSumType*ThirdSumType)/(SecondSumType*SecondSumType);
    }
    }
 all = FirstTerm - SecondTerm;
    //all[iEta] = FirstTerm - SecondTerm;
//  }
    return all;
  }



