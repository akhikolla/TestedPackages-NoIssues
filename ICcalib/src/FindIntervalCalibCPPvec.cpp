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
NumericVector FindIntervalCalibCPPvec(NumericVector w, NumericVector wres) {
  int npoints = w.size();
  NumericVector out(2);
  bool CondLocation;
      int j=0;
      CondLocation = true;
      while(CondLocation == true)
      {
      if (wres[j] == 1)
      {
        if(j==0)
        {
          out[0] = 0;
          out[1] = w(0);
          CondLocation = false;
        } else
          {     
          out[0] = w[j-1];
          out[1] = w[j];
          CondLocation = false;
          }
      } else if (wres[j]==INFINITY)
      {
        if (j==0)
        {
          out[0] = 0;
          out[1] = INFINITY;
        } else {
        out[0] = w[j-1];
        out[1] = INFINITY;
        }
        CondLocation = false;
      } else if (j==npoints-1) {
        out[0] = w[j];
        out[1] = INFINITY;
        CondLocation = false;}
      j += 1;
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
