//Includes/namespaces
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/distributions/normal.hpp>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

boost::math::normal dist(0.0, 1.0);

// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) 
{
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
void normalizeInput(NumericVector& x) 
{
  IntegerVector rank = match(x, stl_sort(x));
  int m = x.size()+1;
  for(int i = 0; i < x.size(); ++i)
  {
    double p = (double)rank[i]/m;
    x[i] = quantile(dist, p);
  }
}

//' Normalize the input Matrix
//' 
//' @param x A numeric matrix
//' @export
// [[Rcpp::export]]
NumericMatrix normalizeInputMatrix(NumericMatrix& x)
{
  int nrow = x.nrow();
  for (int i = 0; i < nrow; i++) {
    NumericVector v =  x.row(i);
    normalizeInput(v);
    x.row(i) = v;
  }
  return x;
}






