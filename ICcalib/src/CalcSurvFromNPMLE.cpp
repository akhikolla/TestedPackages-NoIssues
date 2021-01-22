#include <Rcpp.h>
using namespace Rcpp;


// Given a 2-row matrix, where each column is the turnbull interval, and a vecotr of probabilities
// return the estimated probability 

//data is assumed to be sorted by time
// [[Rcpp::export]]
NumericVector CalcSurvFromNPMLE(NumericVector probs, NumericVector points, NumericMatrix Tbull) {
  int ncol = Tbull.ncol();
  int nSurvProbs = points.size();
  int j=0;
  NumericVector SurvProbs(nSurvProbs);

  bool CondLocation;
  for (int i = 0; i < nSurvProbs; ++i)
  {
    CondLocation=true;
    SurvProbs[i] = 1;
    if (points[i] > Tbull(0,0))
    {
      j=0;
      while(CondLocation == true)
      {
        if (points[i] < Tbull(1,j))
        {
          if (points[i] > Tbull(0,j))
          {
            SurvProbs[i] -= probs[j]*(points[i]-Tbull(0,j))/(Tbull(1,j)-Tbull(0,j));
            CondLocation = false;
          }
          else 
          {
            CondLocation = false;
          }
        }
        else {
          SurvProbs[i] -= probs[j];
          if (j==ncol-1)
          {
            CondLocation = false;
          }
          j +=1;
        }
      }
    }
  }
  return SurvProbs;
}


