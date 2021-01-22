#ifndef hpa_hpaSelection_H
#define hpa_hpaSelection_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

Rcpp::List hpaSelection(Rcpp::Formula selection,
	Rcpp::Formula outcome,
	Rcpp::DataFrame data,
	int z_K,
	int y_K,
	int pol_elements,
	bool is_Newey,
	Rcpp::NumericVector x0,
	bool is_Newey_loocv,
	double z_sd_fixed,
	bool is_parallel,
  bool is_validation);

Rcpp::List hpaSelectionLnLOptim_List(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

double hpaSelectionLnLOptim(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

Rcpp::NumericVector hpaSelectionLnLOptim_ind(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

List hpaSelectionLnLOptim_grad_List(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

NumericVector hpaSelectionLnLOptim_grad(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

NumericMatrix hpaSelectionLnLOptim_grad_ind(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

NumericMatrix hpaSelectionLnLOptim_hessian(Rcpp::NumericVector x0, Rcpp::List hpaSelection_args);

Rcpp::List predict_hpaSelection(Rcpp::List object, 
	Rcpp::DataFrame newdata, 
	std::string method, 
	bool is_cond,
	bool is_outcome);

double logLik_hpaSelection(Rcpp::List model);

Rcpp::List plot_hpaSelection(Rcpp::List x, bool is_outcome);

void print_summary_hpaSelection(Rcpp::List x);

Rcpp::List summary_hpaSelection(Rcpp::List object);

#endif
