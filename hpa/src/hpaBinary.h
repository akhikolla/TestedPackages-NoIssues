#ifndef hpa_hpaBinary_H
#define hpa_hpaBinary_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

Rcpp::List hpaBinary(Rcpp::Formula formula,
	Rcpp::DataFrame data,
	int K,
	double z_mean_fixed,
	double z_sd_fixed,
	double z_constant_fixed,
	bool is_z_coef_first_fixed,
	bool is_x0_probit,
	bool is_sequence,
	Rcpp::NumericVector x0,
	Rcpp::String cov_type,
	int boot_iter,
	bool is_parallel,
	String opt_type,
	Rcpp::List opt_control,
	bool is_validation);

List hpaBinaryLnLOptim_List(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

double hpaBinaryLnLOptim(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

NumericVector hpaBinaryLnLOptim_ind(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

Rcpp::List hpaBinaryLnLOptim_grad_List(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

Rcpp::NumericVector hpaBinaryLnLOptim_grad(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

Rcpp::NumericMatrix hpaBinaryLnLOptim_grad_ind(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

Rcpp::NumericMatrix hpaBinaryLnLOptim_hessian(Rcpp::NumericVector x0, Rcpp::List hpaBinary_args);

Rcpp::NumericVector predict_hpaBinary(Rcpp::List object, 
	Rcpp::DataFrame newdata, 
	bool is_prob);

Rcpp::List summary_hpaBinary(Rcpp::List object);

void print_summary_hpaBinary(Rcpp::List x);

void plot_hpaBinary(Rcpp::List x);

double logLik_hpaBinary(Rcpp::List object);

#endif
