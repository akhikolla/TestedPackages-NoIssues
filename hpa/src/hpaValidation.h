#ifndef hpa_hpaValidation_H
#define hpa_hpaValidation_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

void pol_Validate(NumericVector pol_degrees,
                  NumericVector pol_coefficients);

void ind_Validate(LogicalVector given_ind,
                  LogicalVector omit_ind);

void mean_Validate(NumericVector mean);
void sd_Validate(NumericVector sd);

void expectation_powers_Validate(NumericVector sd);

#endif
