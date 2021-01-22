#include <Rcpp.h>
#include "Emcdf.h"

using namespace Rcpp;

double logR(const NumericMatrix& table, const int& n, const int& i, const int& j){

  int n1 = table.at(j, i);
  int n2 = table.at(n-1, i) - n1;
  int n3 = table.at(j, n-1) - n1;
  int n4 = n - n1 - n2 - n3;

  int ox = table.at(n-1, i);
  int oy = table.at(j, n-1);
  int rx = n - ox;
  int ry = n - oy;

  double term1 = 0;
  double term2 = 0;
  double term3 = 0;
  double term4 = 0;


  if(n1 != 0)
    term1 = n1*log((double)(ox*oy)/(n1*n));
  else
    term1 = 0;
  if(n2 != 0)
    term2 = n2*log((double)(ox*ry)/(n2*n));
  else
    term2 = 0;
  if(n3 != 0)
    term3 = n3*log((double)(rx*oy)/(n3*n));
  else
    term3 = 0;
  if(n4 != 0)
    term4 = n4*log((double)(rx*ry)/(n4*n));
  else
    term4 = 0;


   return term1 + term2 + term3 + term4;
}

// [[Rcpp::export]]
double Tn(NumericVector& x, NumericVector& y){
  int n = x.length();
  Emcdf emcdf(x, y, false);
  NumericMatrix& table = emcdf.getTable();

  double result = 0;

  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      result += logR(table, n, i, j);

  return (-2*result)/(n*n);

}

// [[Rcpp::export]]
int MC_EL_count(double tn, int n, int sn){
  int count = 0;
  for(int i=0; i<sn; ++i){
    NumericVector a = rnorm(n);
    NumericVector b = rnorm(n);
    if(tn < Tn(a, b))
      ++count;
  }
  return count;
}

// [[Rcpp::export]]
NumericVector randEl(int& n, int& sn){
  std::vector<double> out;
  for(int i=0; i<sn; ++i){
    NumericVector x = rnorm(n);
    NumericVector y = rnorm(n);
    out.push_back(Tn(x, y));
  }
  return wrap(out);
}
