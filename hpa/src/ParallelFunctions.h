#ifndef hpa_ParallelFunctions_H
#define hpa_ParallelFunctions_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

NumericVector ParallelVectorPow(NumericVector x, double value);

NumericVector ParallelVectorExp(NumericVector x);

NumericVector dnorm_parallel (NumericVector x, double mean, 
							  double sd, bool is_parallel);

NumericVector pnorm_parallel(NumericVector x, double mean, 
							 double sd, bool is_parallel);

#endif
