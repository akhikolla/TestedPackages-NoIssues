#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector FindIntervalCPP(double point, NumericMatrix w) {
  int nrow = w.nrow();
  int ncol = w.ncol();
  
  IntegerVector intervalW(nrow);
  bool CondLocation;
    for (int i = 0; i < nrow; ++i)
  {
      int j=0;
      CondLocation = true;
      while(CondLocation == true)
      {
      if (w(i,j) > point)
      {
        if(j==0)
        {
          intervalW[i] = 1;
          CondLocation = false;
        } else
          {     
          intervalW[i] = j + 1;
          CondLocation = false;
          }
      } else if (j==ncol-1) {
        intervalW[i] = ncol+1;
        CondLocation = false;}
      j += 1;
      }
  }
        
    return intervalW;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
