#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List reg(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	NumericMatrix x, 
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericMatrix beta,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, const int m,
	double betasl, NumericVector xplot1,NumericVector xplot2,
  	NumericVector time, int thin);

// [[Rcpp::export]]
List reg_resume(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	NumericMatrix x, 
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericMatrix beta,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, const int m,
	double betasl, NumericVector xplot1,NumericVector xplot2,
  	NumericVector time, int thin, std::vector<int> emptybasket, IntegerVector allbaskets);

// [[Rcpp::export]]
List predreg(NumericMatrix alpharec,
       NumericMatrix lambdarec,
       NumericMatrix betarec,
       NumericMatrix xplot, NumericVector tplot, double alpha);
