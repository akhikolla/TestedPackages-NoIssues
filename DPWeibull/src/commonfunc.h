#include<Rcpp.h>
using namespace Rcpp ;
//find the percentile of each column in a matrix
NumericVector colpercentileRcpp(NumericMatrix x, double percentile);

NumericVector matrixtimesvector(NumericMatrix x, NumericVector beta);

//make a submatrix by making the unwanted column be 0 
NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition);
//sample from a multinomial distribution
IntegerVector oneMultinomCalt(NumericVector probs) ;

//concentration parameter generation, use the method from Escobar and West
double nugen(const double nu, const int npts, const int ngrp,const double a, const double b);

// for fixed lambda, find the lower boundary for alpha
double findbase(double lambda);
// for fixed alpha, find the lower boundary for log lambda
double inversebase(double alpha);
//if the result is a real value
bool testreal(double d);


double truncauchy(double truncauchy, void* truncauchy_data);
void samptruncauchy(double *xsamp,double sl);

