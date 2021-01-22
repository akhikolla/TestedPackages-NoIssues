
#include <Rcpp.h>
#include <R.h>
#include <algorithm>
#include <vector>
#include "Emcdf.h"

using namespace Rcpp;

double EmF(const NumericMatrix& table, int n, int i, int j){
  if(i>n-1) i = n-1;
  if(j>n-1) j = n-1;
  if(i<=0 && j<=0){return 0.25*table.at(0,0)/n;}
  if(i<=0){return (0.5*table.at(j-1,0) + 0.25*(table.at(j, 0) - table.at(j-1, 0)))/n;}
  if(j<=0){return (0.5*table.at(0, i-1) + 0.25*(table.at(0, i) - table.at(0, i-1)))/n;}

  double pp = table.at(j, i);
  double pm = table.at(j, i-1);
  double mp = table.at(j-1, i);
  double mm = table.at(j-1, i-1);

  return (mm + 0.5*(mp - mm) + 0.5*(pm - mm) + 0.25*(mm+pp-mp-pm))/n;

}

NumericVector* rank(const NumericVector& x){
  int n = x.length();
  NumericVector* out = new NumericVector(n);
  double cur = 0;
  for(int i=0; i<n; ++i){
    cur = x.at(i);
    for(int j=0; j<n; ++j){
      if(cur >= x.at(j))
        ++(out->at(i));
    }
  }
  return out;
}


// [[Rcpp::export]]
double vex(NumericVector& x, NumericVector& y){

  int n = x.length();
  Emcdf emcdf(x,y, false);
  //x sorted, y sorted by x,

  NumericMatrix table = emcdf.getTable();
  const int r = round(0.5*pow(n,0.8));
  const int& m = r;

  double delta = 0;

  NumericVector* rankY = rank(y);
  int I = 0;
  for(int Si=0; Si<n; ++Si){
    I = rankY->at(Si) - 1;
    double F1 = EmF(table, n, Si + r, I + m);
    double F2 = EmF(table, n, Si - r, I + m);
    double F3 = EmF(table, n, Si + r, I - m);
    double F4 = EmF(table, n, Si - r, I - m);

    double F5 = 0;
    double F6 = 0;

    if(Si + r < n)
      F5 = (double)(Si + r + 1)/n;
    else
      F5 = 1;

    if(Si - r >= 0)
      F6 = (double)(Si - r + 1)/n;
    else
      F6 = (double)1/n;

    double den = 0;
    if(F5 - F6 ==0){den = (double)1/n;}
    else{den = F5 - F6;}

    delta += std::log((F1-F2-F3+F4+pow(n,-0.45))/den);
  }

  delete rankY;

  return 0.2*n*std::log(n) + delta;

}

// [[Rcpp::export]]
int MC_count(double Ts, int n, int sn){
  int count = 0;
  for(int i=0; i<sn; ++i){
    NumericVector a = rnorm(n);
    NumericVector b = rnorm(n);
    if(vex(a,b) > Ts)
      ++count;
  }

  return count;
}

// [[Rcpp::export]]
NumericVector randTs(int& n, int& sn){
  std::vector<double> out;
  for(int i=0; i<sn; ++i){
    NumericVector x = rnorm(n);
    NumericVector y = rnorm(n);
    out.push_back(vex(x, y));
  }
  return wrap(out);
}



