#ifndef hpa_hpaML_H
#define hpa_hpaML_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

List hpaML(NumericMatrix x,
	NumericVector pol_degrees,
	NumericVector tr_left,
	NumericVector tr_right,
	LogicalVector given_ind,
	LogicalVector omit_ind,
	NumericVector x0,
	String cov_type,
	int bootstrap,
	bool is_parallel,
	String opt_type,
	List opt_control,
	bool is_validation);

List hpaLnLOptim_List(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

double hpaLnLOptim(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericVector hpaLnLOptim_ind(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

List hpaLnLOptim_grad_List(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericVector hpaLnLOptim_grad(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericMatrix hpaLnLOptim_grad_ind(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

NumericMatrix hpaLnLOptim_hessian(Rcpp::NumericVector x0, Rcpp::List hpaML_args);

Rcpp::NumericVector predict_hpaML(Rcpp::List object,
	Rcpp::NumericMatrix newdata);

Rcpp::NumericVector mecdf(NumericMatrix x);

void print_summary_hpaML(Rcpp::List x);

Rcpp::List summary_hpaML(Rcpp::List model);

double logLik_hpaML(Rcpp::List model);

Rcpp::StringVector starVector(Rcpp::NumericVector p_values);

#endif
