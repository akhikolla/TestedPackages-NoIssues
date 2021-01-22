/*****************************/
/* Jiaxing Lin               */
/* 05-13-2016                */
/* Estimate Alpha without    */
/*****************************/

#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::NumericVector alphaEst1(
    Rcpp::NumericVector dtimeFactor,
    Rcpp::NumericVector dtime,
    Rcpp::NumericVector delta) 
{
  double d = 0;
  double r = 0;
  
  int NPat = dtime.size();
  
  int m = dtimeFactor.size();
  Rcpp::NumericVector alpha(m);
  
  for(int j=0; j < m; j++){
    d = 0;
    r = 0;
    for(int i=0; i < NPat; i++)
    {
      if(dtime[i] == dtimeFactor[j] && delta[i] == 1)
        d++;
      if(dtime[i] > dtimeFactor[j])
        r++; 	
    } 
    if(r+d != 0)	
      alpha[j] = r/(r+d);
    else
      alpha[j] = 0;
  }
  
  return alpha; 
}




