
#include <Rcpp.h>
#include <R.h>
#include <vector>
#include "Emcdf.h"

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix emcdf_output(NumericVector& x, NumericVector& y, bool is_tie){

  Emcdf emcdf(x,y, is_tie);
  return(emcdf.getTable());

}
