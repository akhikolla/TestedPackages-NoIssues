#include <Rcpp.h>
#include <complex>
#include <cmath>
#include "QF.h"

using namespace Rcpp;



// [[Rcpp::export]]
std::vector<double> dQF_c(std::vector<double> q, List Mellin_list){
  // take elements from the Mellin
  std::vector<std::complex<double> > Mellin= Rcpp::as<std::vector<std::complex<double> > >( Mellin_list["Mellin"]);
  std::vector<std::complex<double> > z= Rcpp::as<std::vector<std::complex<double> > >(Mellin_list["z"]);
  std::complex<double> delta = Rcpp::as<std::complex<double> >(Mellin_list["delta"]);
  // define quantities
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> f_q_prov, pii(0.0, 2.0 * M_PI);
  std::vector<double> f_q(dim_q);
  // compute densities
  for(int j = 0; j < dim_q; j++){
    f_q_prov = 0;
    for(int k = 0; k < dim_Mellin; k++){
      f_q_prov += Mellin[k] * pow(q[j], - z[k]) * delta;
    }
    f_q_prov = f_q_prov / pii;
    f_q[j] = f_q_prov.real();
  }

  return(f_q);

}

// [[Rcpp::export]]
std::vector<double> dQF_c_scal(std::vector<double> q, List Mellin_list){
  // take elements from the Mellin
  std::vector<std::complex<double> > Mellin= Rcpp::as<std::vector<std::complex<double> > >( Mellin_list["Mellin"]);
  std::vector<std::complex<double> > z= Rcpp::as<std::vector<std::complex<double> > >(Mellin_list["z"]);
  std::complex<double> delta = Rcpp::as<std::complex<double> >(Mellin_list["delta"]);
  double lambda_min = Mellin_list["lambda_min"];
  // define quantities
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> f_q_prov, pii(0.0, 2.0 * M_PI);
  std::vector<double> f_q(dim_q);
  // compute densities
  for(int j = 0; j < dim_q; j++){
    f_q_prov = 0;
    for(int k = 0; k < dim_Mellin; k++){
      f_q_prov += Mellin[k] * pow(q[j] / lambda_min, - z[k]) * delta;
    }
    f_q_prov = f_q_prov / pii;
    f_q[j] = f_q_prov.real() / lambda_min;
  }

  return(f_q);

}


// [[Rcpp::export]]
std::vector<double> pQF_c(std::vector<double> q, List Mellin_list){
  // take elements from the Mellin
  std::vector<std::complex<double> > Mellin = Rcpp::as<std::vector<std::complex<double> > >( Mellin_list["Mellin"]);
  std::vector<std::complex<double> > z = Rcpp::as<std::vector<std::complex<double> > >(Mellin_list["z"]);
  std::complex<double> delta = Rcpp::as<std::complex<double> >(Mellin_list["delta"]);
  // define quantities
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> F_q_prov, pii(0.0, 2.0 * M_PI);
  std::vector<double> F_q(dim_q);
  // compute cumulative functions
  for(int j = 0; j < dim_q; j++){
    F_q_prov = 0.0;
    for(int k = 0; k < dim_Mellin; k++){
      F_q_prov += - Mellin[k] / (z[k] - 1.0) * pow(q[j], - z[k] + 1.0) * delta;
    }
    F_q_prov = F_q_prov / pii;
    F_q[j] = F_q_prov.real();
  }
  return(F_q);
}


// [[Rcpp::export]]
std::vector<double> pQF_c_scal(std::vector<double> q, List Mellin_list){
  // take elements from the Mellin
  std::vector<std::complex<double> > Mellin = Rcpp::as<std::vector<std::complex<double> > >( Mellin_list["Mellin"]);
  std::vector<std::complex<double> > z = Rcpp::as<std::vector<std::complex<double> > >(Mellin_list["z"]);
  std::complex<double> delta = Rcpp::as<std::complex<double> >(Mellin_list["delta"]);
  double lambda_min = Mellin_list["lambda_min"];
  // define quantities
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> F_q_prov, pii(0.0, 2.0 * M_PI);
  std::vector<double> F_q(dim_q);
  // compute cumulative functions
  for(int j = 0; j < dim_q; j++){
    F_q_prov = 0.0;
    for(int k = 0; k < dim_Mellin; k++){
      F_q_prov += - Mellin[k] / (z[k] - 1.0) * pow(q[j] / lambda_min, - z[k] + 1.0) * delta;
    }
    F_q_prov = F_q_prov / pii;
    F_q[j] = F_q_prov.real();
  }
  return(F_q);
}
