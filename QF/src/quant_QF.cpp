#include <RcppGSL.h>
#include <complex>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> dQF_quant(std::vector<double> q, std::vector<std::complex<double> > Mellin,
                              std::vector<std::complex<double> > z,
                              std::complex<double> delta){
  // define quantites
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> f_q_comp, pii(0.0, 2.0 * M_PI);
  std::vector<double> f_q(dim_q);
  // compute density
  for(int j=0; j < dim_q; j++){
    f_q_comp = 0.0;
    for(int k = 0; k < dim_Mellin; k++){
      f_q_comp += Mellin[k] * pow(q[j], - z[k]) * delta;
    }
    f_q_comp = f_q_comp / pii;
    f_q[j] = f_q_comp.real();
  }
  return(f_q);
}

// [[Rcpp::export]]
std::vector<double> pQF_quant(std::vector<double> q, std::vector<std::complex<double> > Mellin,
                              std::vector<std::complex<double> > z,
                              std::complex<double> delta){
  // define quantites
  int dim_q = q.size(), dim_Mellin = Mellin.size();
  std::complex<double> F_q_comp, pii(0.0, 2.0 * M_PI);
  std::vector<double> F_q(dim_q);
  // compute cumulative distribution
  for(int j=0; j < dim_q; j++){
    F_q_comp = 0.0;
    for(int k = 0; k < dim_Mellin; k++){
      F_q_comp += - Mellin[k] / (z[k] - 1.0)* pow(q[j], - z[k] + 1.0) * delta;
    }
    F_q_comp = F_q_comp / pii;
    F_q[j] = F_q_comp.real();
  }
  return(F_q);
}

// [[Rcpp::export]]
std::vector<double> qQF_c(std::vector<double> p, List Mellin_list,
                          double eps_quant, int maxit_quant, double q0) {
  // define quantities
  int n = p.size();
  std::vector<double> q_p(n), q_p_new(1), q_p_old(1), F_q(maxit_quant), f_q(maxit_quant);
  double error_quant = 99.0;
  // take elements from the Mellin
  std::vector<std::complex<double> > Mellin = Rcpp::as<std::vector<std::complex<double> > >( Mellin_list["Mellin"]);
  std::vector<std::complex<double> > z = Rcpp::as<std::vector<std::complex<double> > >(Mellin_list["z"]);
  std::complex<double> delta = Rcpp::as<std::complex<double> >(Mellin_list["delta"]);
  // compute quantiles
  for(int i = 0; i < n; i++){
    error_quant = 99.0;
    q_p_new[0] = q0;
    for(int k=0; k<maxit_quant; k++){
      F_q[k] = pQF_quant(q_p_new, Mellin, z, delta)[0];
      f_q[k] = dQF_quant(q_p_new, Mellin, z, delta)[0];
      q_p_old = q_p_new;
      q_p_new[0] = q_p_new[0] - (F_q[k]- p[i]) / f_q[k];
      //Rprintf("f: %lf; F: %Lf \n", f_q[k], F_q[k]);
      while(q_p_new[0] < 0.0 || f_q[k] < 0.0 || F_q[k] > 1.0){
        q_p_old[0] = q_p_old[0] / 2.0;
        F_q[k] = pQF_quant(q_p_old, Mellin, z, delta)[0];
        f_q[k] = dQF_quant(q_p_old, Mellin, z, delta)[0];
        q_p_new[0] = q_p_old[0] - (F_q[k] - p[i]) / f_q[k];
        //Rprintf("f: %lf; F: %Lf \n", f_q[k], F_q[k]);
      }
      error_quant= fabs((F_q[k] - p[i]) / p[i]) ;
      q_p[i] = q_p_new[0];
      if(error_quant < eps_quant) break;
      if(k == maxit_quant-1){
        warning("Expected precision not reached in computing quantiles");
      }

    }
  }
  return q_p;
}

