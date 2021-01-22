#include <RcppArmadillo.h>
using namespace Rcpp;


// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// beta- a value

//data is not assumed to be sorted by time

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double CoxLogLik(arma::vec betagamma, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Z) {
  int n = tm.size();
  int nGamma = betagamma.size()-1;
  double denom=0;
  double logDenom=0;
  double logNumer=0;
  double logLik=0;
  int iCaseNum=-1;
  double beta = betagamma[0];
  arma::vec gamma = betagamma.subvec(1,nGamma);
  arma::mat contrib=1 + ps*(exp(beta)-1);
  arma::vec GamZ = Z * gamma;
  arma::vec ExpGamZ = exp(Z * gamma);
  for (int i = 0; i < n; ++i)
  {
    if (event[i]) {
      iCaseNum += 1;
      logNumer += log(contrib(iCaseNum,i)) + GamZ[i];
      denom = contrib(iCaseNum,i)*ExpGamZ[i];
      for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        denom += contrib(iCaseNum,j)*ExpGamZ[j];
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


