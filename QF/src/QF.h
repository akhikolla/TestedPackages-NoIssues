#ifndef _QF_H
#define _QF_H

#include <Rcpp.h>
#include <complex>
#include <R_ext/Rdynload.h>

using namespace Rcpp;

std::complex<double> pochhammer_complex(std::complex<double> x, std::complex<double> y);

std::vector<double> compute_ak_c(NumericVector lambdas,
                                 int maxit,
                                 double eps,
                                 double beta);
std::vector<double> compute_ak_nc(NumericVector lambdas,
                                  NumericVector etas,
                                  int maxit,
                                  double eps,
                                  double beta);

std::vector<double> dQF_c(std::vector<double> q, List Mellin_list);
std::vector<double> dQF_c_scal(std::vector<double> q, List Mellin_list);
std::vector<double> pQF_c(std::vector<double> q, List Mellin_list);
std::vector<double> pQF_c_scal(std::vector<double> q, List Mellin_list);


List get_mellin_QF(NumericVector lambdas,
                   NumericVector etas,
                   double rho, double h,
                   std::complex<double> delta,
                   double eps, double eps_quant,
                   int maxit_ak, int maxit_quant,
                   int maxit_delta, double step_delta);

List get_mellin_QF_ratio(NumericVector lambdas_1, NumericVector lambdas_2,
                         NumericVector etas_1, NumericVector etas_2,
                         double rho, double h,
                         std::complex<double> delta,
                         double eps, double eps_quant,
                         int maxit_ak, int maxit_quant,
                         int maxit_delta, double step_delta);

List Mellin_QF(NumericVector lambdas,
               std::vector<double> a_k,
               int maxit,
               double h, std::complex<double> delta,
               double eps, double beta);

List Mellin_QF_error(NumericVector lambdas,
                     std::vector<double> a_k,
                     int maxit,
                     double h, std::complex<double> delta,
                     double eps, double beta,
                     std::vector<double> q_lims);

List Mellin_QF_ratio(NumericVector lambdas_1, NumericVector lambdas_2,
                     std::vector<double> a_k_1, std::vector<double> a_k_2,
                     int maxit,
                     double h, std::complex<double> delta,
                     double eps, double beta_1, double beta_2);

List Mellin_QF_ratio_error(NumericVector lambdas_1, NumericVector lambdas_2,
                           std::vector<double> a_k_1, std::vector<double> a_k_2,
                           int maxit,
                           double h, std::complex<double> delta,
                           double eps, double beta_1, double beta_2,
                           std::vector<double> q_lims);

std::vector<double> qQF_c(std::vector<double> p, List Mellin_list,
                          double eps_quant, int maxit_quant, double q0);

double find_maximum_error(std::vector<double> range_f_ref,
                          std::vector<double> range_F_ref,
                          std::vector<double> range_f_new,
                          std::vector<double> range_F_new);

#endif
