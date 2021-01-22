#include <Rcpp.h>
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

//data is assumed to be sorted by time
// [[Rcpp::export]]
NumericMatrix FindIntervalCalibCPP(NumericMatrix w, NumericMatrix wres) {
  int nrow = w.nrow();
  int ncol = w.ncol();
  NumericMatrix out(nrow,2);
  bool CondLocation;
    for (int i = 0; i < nrow; ++i)
  {
      int j=0;
      CondLocation = true;
      while(CondLocation == true)
      {
      if (wres(i,j) == 1)
      {
        if(j==0)
        {
          out(i,0) = 0;
          out(i,1) = w(i,0);
          CondLocation = false;
        } else
          {     
          out(i,0) = w(i,j-1);
          out(i,1) = w(i,j);
          CondLocation = false;
          }
      } else if (wres(i,j)==INFINITY)
      {
        if (j==0)
        {
          out(i,0) = 0;
          out(i,1) = INFINITY;
        } else {
        out(i,0) = w(i,j-1);
        out(i,1) = INFINITY;
        }
        CondLocation = false;
      } else if (j==ncol-1) {
        out(i,0) = w(i,j);
        out(i,1) = INFINITY;
        CondLocation = false;}
      j += 1;
      }
  }
        
    return out;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
