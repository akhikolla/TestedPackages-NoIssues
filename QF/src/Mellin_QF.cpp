#include <RcppGSL.h>
#include <complex>
#include "QF.h"

using namespace Rcpp;


// [[Rcpp::export]]
List Mellin_QF(NumericVector lambdas,
                    std::vector<double> a_k,
                    int maxit,
                    double h, std::complex<double> delta,
                    double eps, double beta) {
  // define quantities
  int r = lambdas.size(), K = a_k.size();
  double alpha = r / 2.0, Tn;
  std::complex<double> Mellin_prov, P_k;
  std::vector<std::complex<double> > Mellin, z, z_compl, Mellin_compl;
  // init vector of evaluation points
  z.push_back(h);
  // Compute Mellin until the desired error is reached
  for (int t = 0; t < maxit; t++) {
    Mellin_prov = 0.0;
    P_k = 1.0;
    for (int k = 0; k < K; k++) {
      Mellin_prov += a_k[k] * P_k;
      P_k *= 1.0 + (z[t] - 1.0) / (k * 1.0 + alpha);
    }
    Mellin_prov = (Mellin_prov / (pow((2.0 * beta), 1.0 - z[t]))) *
      pochhammer_complex(alpha, z[t] - 1.0);
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
List Mellin_QF_error(NumericVector lambdas,
                      std::vector<double> a_k,
                      int maxit,
                      double h, std::complex<double> delta,
                      double eps, double beta,
                      std::vector<double> q_lims) {
  // define quantities
  int r = lambdas.size(), K = a_k.size();
  double alpha = r / 2.0, Tn;
  std::complex<double> Mellin_prov, P_k;
  std::vector<std::complex<double> > Mellin, z, z_compl, Mellin_compl;
  // init vector of evaluation points
  z.push_back(h);
  //Rprintf("low %lf, up %lf\n", q_lims[0], q_lims[1]);

  // Compute Mellin until the desired error is reached
  for (int t = 0; t < maxit; t++) {
    Mellin_prov = 0.0;
    P_k = 1.0;
    for (int k = 0; k < K; k++) {
      Mellin_prov += a_k[k] * P_k;
      P_k *= 1.0 + (z[t] - 1.0) / (k * 1.0 + alpha);
    }
    Mellin_prov = (Mellin_prov / (pow((2.0 * beta), 1.0 - z[t]))) *
      pochhammer_complex(alpha, z[t] - 1.0);
    Mellin.push_back(Mellin_prov);
    Tn = z[t].imag();
    z.push_back(z[t] + delta);
    // Check error
    //Rprintf("dens: %lf\n rip: %lf\n", std::abs(Mellin_prov) * (h * h + Tn * Tn) * (M_PI_2 - atan(Tn / h)) / (M_PI * h * pow(q_lims[0], h)),
    //       std::abs(Mellin_prov) * pow(q_lims[1], 1.0-h) * (h * h + Tn * Tn) * (1.0 - Tn / pow(Tn * Tn+ pow(h-1.0, 2.0),0.5)) / (M_PI * pow(h-1.0,2.0)));

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


// [[Rcpp::export]]
double find_maximum_error(std::vector<double> range_f_ref,
                          std::vector<double> range_F_ref,
                          std::vector<double> range_f_new,
                          std::vector<double> range_F_new) {
  std::vector<double> err_extr_f(2), err_extr_F(2);
  double err_max_f, err_max_F, err_max;

  for(int k=0; k<2; k++){
    err_extr_f[k] = std::fabs(range_f_ref[k]-range_f_new[k]);///range_f_ref[k];
    err_extr_F[k] = std::fabs(range_F_ref[k]-range_F_new[k]);///range_F_ref[k];
  }
  if(err_extr_f[0]<err_extr_f[1]){
    err_max_f = err_extr_f[1];
  }else{
    err_max_f = err_extr_f[0];
  }
  if(err_extr_F[0]<err_extr_F[1]){
    err_max_F = err_extr_F[1];
  }else{
    err_max_F = err_extr_F[0];
  }

  if(err_max_f<err_max_F){
    err_max = err_max_F;
  }else{
    err_max = err_max_f;
  }
  return(err_max);
}
