#include <Rcpp.h>
using namespace Rcpp;


// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// beta- a value

//data is not assumed to be sorted by time
// [[Rcpp::export]]
NumericVector Calcb(double beta, NumericVector tm, LogicalVector event, NumericMatrix ps) {
  int n = tm.size();
  int sumD = sum(event);
  NumericVector b(n);
  NumericVector Szero(sumD);
  NumericVector Sone(sumD);
  int iCaseNum= -1;
  int kCaseNum= -1;
  int mCaseNum= -1;
  NumericMatrix contribPbeta= ps*exp(beta);
//  NumericMatrix contribPbetaMinusP= (1-ps)*ps*exp(beta); Not working!!
  NumericMatrix contribDenom= 1 + ps*(exp(beta)-1);
// First a loop to calculate s0 and s1
    for (int i = 0; i < n; ++i)
  {
    if (event[i]) {
      iCaseNum += 1;
      for(int j = 0; j < n; ++j) {
        if (tm[j]>tm[i]) {
      Szero(iCaseNum) += contribDenom(iCaseNum,j);
      Sone(iCaseNum) += contribPbeta(iCaseNum,j);
        }
      }
    }}
    for (int k = 0; k < n; ++k)
    {
      if (event[k]) {
        kCaseNum += 1;
      b(k) += contribPbeta(kCaseNum,k)/contribDenom(kCaseNum,k) - Sone(kCaseNum)/Szero(kCaseNum);
      }
      mCaseNum = -1;
      for(int m = 0; m < n; ++m) {
     if (event[m]) {
         mCaseNum +=1;
        if (tm[k]>tm[m]) {
          b(k) -= (n*contribDenom(mCaseNum,k)/Szero(mCaseNum)) * (contribPbeta(mCaseNum,k)/contribDenom(mCaseNum,k) - 
            Sone(mCaseNum)/Szero(mCaseNum))/n;
    }}}}   
    return b;
  }

