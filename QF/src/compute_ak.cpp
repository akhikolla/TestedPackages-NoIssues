#include <Rcpp.h>
#include "QF.h"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> compute_ak_nc(NumericVector lambdas,
                               NumericVector etas,
                               int maxit,
                               double eps,
                               double beta){
  // Inits quantities
  int r = lambdas.size();
  double check_error, eps_a_0, a_0 = 0.0, b_k_par1, b_k_par2;
  double eta = Rcpp::sum(etas);
  std::vector<double> b_k, a_k, b_k_const(r), c_i(r);
  a_k.reserve(maxit);
  b_k.reserve(maxit - 1);

  // create useful quantites
  for (int i = 0; i < r; i++){
    a_0 += std::log(lambdas[i]);
    c_i[i] = 1.0 - beta / lambdas[i];
    b_k_const[i] = etas[i] / lambdas[i];
  };
  // first coefficient
  a_0 = std::exp(0.5 * (- eta + r * std::log(beta) - a_0));
  a_k.push_back(1.0);

  // Quantities for error control
  check_error = 1.0 / a_0 - 1.0;
  eps_a_0 = eps / a_0;

  // loop for the ak computation
  for(int k = 1; k < maxit; k++) {
    b_k_par1 = 0.0;
    b_k_par2 = 0.0;
    //loop for b_k
    for (int j = 0; j < r; j++) {
      b_k_par1 += b_k_const[j] * std::pow(c_i[j], 1.0 * k - 1.0);
      b_k_par2 += pow(c_i[j], 1.0 * k);
    }
    b_k.push_back((1.0 * k) * beta * b_k_par1 + b_k_par2);
    //a_k
    a_k.push_back(0.0);
    for (int l = 0;l < k; l++){
      a_k[k] += b_k[k - l - 1] * a_k[l];
    }
    a_k[k] = a_k[k] / (2.0 * k);
    // check error
    check_error -= a_k[k];
    if (std::fabs(check_error) < eps_a_0) {
      break;
      }
  }
  int K = a_k.size();
  for(int i = 0; i < K; i++) {
    a_k[i] = a_k[i] * a_0;
    }

  return a_k;
}



// [[Rcpp::export]]
std::vector<double> compute_ak_c(NumericVector lambdas,
                                  int maxit,
                                  double eps,
                                  double beta){
  // Inits quantities
  int r = lambdas.size();
  double check_error, eps_a_0, a_0 = 0.0, b_k_par;
  std::vector<double> b_k, a_k, c_i(r);
  a_k.reserve(maxit);
  b_k.reserve(maxit - 1);

  // create useful quantites
  for (int i = 0; i < r; i++){
    a_0 += std::log(lambdas[i]);
    c_i[i] = 1.0 - beta / lambdas[i];
    };
  // first coefficient
  a_0 = exp(0.5 * (r * std::log(beta) - a_0));
  a_k.push_back(1.0);

  // Quantities for error control
  check_error = 1.0 / a_0 - 1.0;
  eps_a_0 = eps / a_0;

  // loop for the ak computation
  for(int k = 1; k < maxit; k++) {
    b_k_par = 0.0;
    //loop for b_k
    for (int j = 0; j < r; j++) {
      b_k_par += std::pow(c_i[j], 1.0 * k);
    }
    b_k.push_back(b_k_par);
    //a_k
    a_k.push_back(0.0);
    for (int l = 0;l < k; l++){
      a_k[k] += b_k[k - l - 1] * a_k[l];
    }
    a_k[k] = a_k[k] / (2.0 * k);
    // check error
    check_error -= a_k[k];
    if (std::fabs(check_error) < eps_a_0) {
      break;
    }
  }
  int K = a_k.size();
  for(int i = 0; i < K; i++) {
    a_k[i] = a_k[i] * a_0;
  }

  return a_k;
}
