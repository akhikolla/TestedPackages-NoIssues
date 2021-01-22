#include <RcppGSL.h>
#include <complex>
#include "QF.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pQF_depratio_c(NumericVector lambdas_1, NumericVector lambdas_2,
                      double h,
                      std::complex<double> delta,
                      double eps,
                      int maxit_ak,
                      int maxit_delta) {

  std::vector<double> ak_1, ak_2, one(1);
  ak_1.reserve(maxit_ak), ak_2.reserve(maxit_ak);
  one[0] = 1.0;
  double F_ref, F_new, err_max;
  double eps_mell_1 = eps *  0.01;
  double eps_mell_2 = eps * 0.01;
  double eps_trunc = eps;
  double eps_discr = eps;


  int length_1=lambdas_1.size();
  int length_2=lambdas_2.size();
  if(length_1==0){
    F_new = 1.0;
  }else if(length_2==0){
    F_new = 0.0;
  }else{
    // compute ak
    ak_1=compute_ak_c(lambdas_1,
                          maxit_ak,
                          eps_mell_1,
                          min(lambdas_1));
    ak_2=compute_ak_c(lambdas_2,
                          maxit_ak,
                          eps_mell_2,
                          min(lambdas_2));

    if(ak_1.size() == maxit_ak || ak_2.size() == maxit_ak){// stop with error if the ak did not converged
      stop("Computation of the a_k coefficients did not converged: consider to increase 'maxit_comp'");
    }


    // initial Mellin_ref
    List Mellin_ref = Mellin_QF_ratio(
      lambdas_1, lambdas_2,  ak_1,  ak_2,
      maxit_ak,
      h, delta,
      eps_trunc, Rcpp::min(lambdas_1), Rcpp::min(lambdas_2));

    F_ref = pQF_c(one, Mellin_ref)[0];


    // first iteration
    delta = std::complex<double> (0.0,0.05) + delta;
    int direction=1;
    List Mellin_new = Mellin_QF_ratio(
      lambdas_1, lambdas_2,  ak_1,  ak_2,
      maxit_ak,
      h, delta,
      eps_trunc, Rcpp::min(lambdas_1), Rcpp::min(lambdas_2));
    F_new = pQF_c(one, Mellin_ref)[0];

    //find maximum error
    err_max = std::fabs(F_new-F_ref);

    // decide direction: 0 if the initial value is too high
    if (err_max > eps_discr){
      direction = 0;
    }

    if(direction == 0){// reduce delta
      delta = delta - std::complex<double>(0.0,0.05);
      for(int i=0; i<maxit_delta; i++){
        delta = delta / 2.0;
        // Compute new quantities
        List Mellin_new = Mellin_QF_ratio(
          lambdas_1, lambdas_2,  ak_1,  ak_2,
          maxit_ak,
          h, delta,
          eps_trunc, Rcpp::min(lambdas_1), Rcpp::min(lambdas_2));
        F_new = pQF_c(one, Mellin_ref)[0];


        //find maximum error
        err_max = std::fabs(F_new-F_ref);

        if (err_max < eps_discr){
          break; //stop if the error is lower that the threshold
        }
        if(i == maxit_delta - 1){// stop with error if maxit_delta is reached without the desired error
          stop("Computation of the integration step 'delta' did not converged: consider to increase 'maxit_comp'");
        }
        // new become reference
        F_ref = F_new;
      }
    }
  }
return F_new;
}

