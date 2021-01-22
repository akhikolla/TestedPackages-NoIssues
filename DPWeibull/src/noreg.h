#include <Rcpp.h>
using namespace Rcpp;
//calculating the log likelihood
double noreg_loglikelihood(const double tl, const double tr,
const int delta, const int pi,
const double lambda, const double alpha);
// Using Neal 8 algorithm to assign cluster to a given observation。
int noreg_group_assign(const double tl, const double tr,
	const int delta,const int pi,
	const int c, const double nu,
	IntegerVector nm,
	NumericVector alpha, NumericVector lambda,
	const double lambda00, const double alpha00, const double alpha0,
	const double alphaalpha, const double alphalambda, const int m,
	IntegerVector allbaskets, std::vector<int> & emptybasket);

void noreg_update(NumericVector tl, NumericVector tr,
	IntegerVector delta, IntegerVector pi,
	IntegerVector c, IntegerVector nm,
	NumericVector alpha, NumericVector lambda, NumericVector lambda0, 
	const double alpha00,const double alpha0,const double lambda00,
	const double alphaalpha, const double alphalambda,
	double* nextngrp);

// [[Rcpp::export]]
List noreg(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector t,
	int m, int thin);


// [[Rcpp::export]]
List noreg_resume(const int burnin, const int iteration,
	NumericVector tl, NumericVector tr,
	IntegerVector delta,
	IntegerVector pi,
	IntegerVector c,
   	IntegerVector nm, 
	NumericVector alpha,
	NumericVector lambda,
	NumericVector lambda0,
	const double alpha00,
	const double alpha0,
	const double lambda00,
	const double alphaalpha,
	const double alphalambda,
	NumericVector nu,
	NumericVector ngrp,
	const double a, const double b,
	const double ymax, NumericVector t,
	int m, int thin, std::vector<int> emptybasket, IntegerVector allbaskets);
