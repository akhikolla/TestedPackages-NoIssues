#include <RcppGSL.h>
#include <complex>
#include "QF.h"
using namespace Rcpp;


// [[Rcpp::export]]
List Mellin_QF_ratio(NumericVector lambdas_1, NumericVector lambdas_2,
                    std::vector<double> a_k_1, std::vector<double> a_k_2,
                    int maxit,
                    double h, std::complex<double> delta,
                    double eps, double beta_1, double beta_2) {
  // define quantities
  int r_1 = lambdas_1.size(), K_1 = a_k_1.size(),
      r_2 = lambdas_2.size(), K_2 = a_k_2.size();
  double alpha_1 = r_1 / 2.0, alpha_2 = r_2 / 2.0, Tn;
  std::complex<double> Mellin_prov, Mellin_prov_1, Mellin_prov_2, z2, P_k;
  std::vector<std::complex<double> > Mellin, z, z_compl, Mellin_compl;
  // init vector of evaluation points
  z.push_back(h);
  // Compute Mellin until the desired error is reached
  for (int t = 0; t < maxit; t++) {
    // compute Mellin 1
    Mellin_prov_1 = 0.0;
    P_k = 1.0;
    for (int k = 0; k < K_1; k++) {
      Mellin_prov_1 += a_k_1[k] * P_k;
      P_k *= 1.0 + (z[t] - 1.0) / (k * 1.0 + alpha_1);
    }
    Mellin_prov_1 = (Mellin_prov_1 / (pow((2.0 * beta_1), 1.0 - z[t]))) *
      pochhammer_complex(alpha_1, z[t] - 1.0);

    // compute Mellin 2
    Mellin_prov_2 = 0.0;
    P_k = 1.0;
    z2 = 2.0 - z[t];
    for (int k = 0; k < K_2; k++) {
      Mellin_prov_2 += a_k_2[k] * P_k;
      P_k *= 1.0 + (z2 - 1.0) / (k * 1.0 + alpha_2);
    }
    Mellin_prov_2 = (Mellin_prov_2 / (pow((2.0 * beta_2), 1.0 - z2))) *
      pochhammer_complex(alpha_2, z2 - 1.0);
    Mellin_prov = Mellin_prov_1 * Mellin_prov_2;


    Mellin.push_back(Mellin_prov);
    Tn = z[t].imag();
    z.push_back(z[t] + delta);
    // Check error
    if(std::abs(Mellin_prov) * (h * h + Tn * Tn) * (M_PI_2 - atan(Tn / h)) / (M_PI * h) < eps) {
      break;
    }
    if(t == maxit - 1){
      stop("Expected precision not reached in computing the Mellin, increase 'maxit_comp'");
    }
  }
  int T = Mellin.size();
  // complete the vectors with the complex conjugates
  for (int t = 1; t < T; t++){
    z_compl.push_back(std::conj(z[T - t]));
    Mellin_compl.push_back(std::conj(Mellin[T - t]));
  }
  for (int t = 0; t < T; t++){
    z_compl.push_back(z[t]);
    Mellin_compl.push_back(Mellin[t]);
  }
  // return
  return Rcpp::List::create(
    Rcpp::Named("Mellin") = Mellin_compl,
    Rcpp::Named("z") = z_compl,
    Rcpp::Named("delta") = delta);

}



// [[Rcpp::export]]
List Mellin_QF_ratio_error(NumericVector lambdas_1, NumericVector lambdas_2,
                     std::vector<double> a_k_1, std::vector<double> a_k_2,
                     int maxit,
                     double h, std::complex<double> delta,
                     double eps, double beta_1, double beta_2,
                     std::vector<double> q_lims) {
  // define quantities
  int r_1 = lambdas_1.size(), K_1 = a_k_1.size(),
    r_2 = lambdas_2.size(), K_2 = a_k_2.size();
  double alpha_1 = r_1 / 2.0, alpha_2 = r_2 / 2.0, Tn;
  std::complex<double> Mellin_prov, Mellin_prov_1, Mellin_prov_2, z2, P_k;
  std::vector<std::complex<double> > Mellin, z, z_compl, Mellin_compl;
  // init vector of evaluation points
  z.push_back(h);
  // Compute Mellin until the desired error is reached
  for (int t = 0; t < maxit; t++) {
    // compute Mellin 1
    Mellin_prov_1 = 0.0;
    P_k = 1.0;
    for (int k = 0; k < K_1; k++) {
      Mellin_prov_1 += a_k_1[k] * P_k;
      P_k *= 1.0 + (z[t] - 1.0) / (k * 1.0 + alpha_1);
    }
    Mellin_prov_1 = (Mellin_prov_1 / (pow((2.0 * beta_1), 1.0 - z[t]))) *
      pochhammer_complex(alpha_1, z[t] - 1.0);

    // compute Mellin 2
    Mellin_prov_2 = 0.0;
    P_k = 1.0;
    z2 = 2.0 - z[t];
    for (int k = 0; k < K_2; k++) {
      Mellin_prov_2 += a_k_2[k] * P_k;
      P_k *= 1.0 + (z2 - 1.0) / (k * 1.0 + alpha_2);
    }
    Mellin_prov_2 = (Mellin_prov_2 / (pow((2.0 * beta_2), 1.0 - z2))) *
      pochhammer_complex(alpha_2, z2 - 1.0);
    Mellin_prov = Mellin_prov_1 * Mellin_prov_2;


    Mellin.push_back(Mellin_prov);
    Tn = z[t].imag();
    z.push_back(z[t] + delta);
    // Check error
    if(std::abs(Mellin_prov) * (h * h + Tn * Tn) * (M_PI_2 - atan(Tn / h)) / (M_PI * h * pow(q_lims[0], h)) < eps &&
       std::abs(Mellin_prov) * pow(q_lims[1], 1.0-h) * (h * h + Tn * Tn) * (1.0 - Tn / pow(Tn * Tn+ pow(h-1.0, 2.0),0.5)) / (M_PI * pow(h-1.0,2.0)) < eps) {
      break;
    }
    if(t == maxit - 1){
      stop("Expected precision not reached in computing the Mellin, increase 'maxit_comp'");
    }
  }
  int T = Mellin.size();
  // complete the vectors with the complex conjugates
  for (int t = 1; t < T; t++){
    z_compl.push_back(std::conj(z[T - t]));
    Mellin_compl.push_back(std::conj(Mellin[T - t]));
  }
  for (int t = 0; t < T; t++){
    z_compl.push_back(z[t]);
    Mellin_compl.push_back(Mellin[t]);
  }
  // return
  return Rcpp::List::create(
    Rcpp::Named("Mellin") = Mellin_compl,
    Rcpp::Named("z") = z_compl,
    Rcpp::Named("delta") = delta);

}
