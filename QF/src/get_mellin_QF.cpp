#include <RcppGSL.h>
#include <complex.h>
#include "QF.h"


// [[Rcpp::export]]
List get_mellin_QF(NumericVector lambdas,
                          NumericVector etas,
                          double rho, double h,
                          std::complex<double> delta,
                          double eps, double eps_quant,
                          int maxit_ak, int maxit_quant,
                          int maxit_delta, double step_delta) {

  std::vector<double> ak,
  range_p(2), range_q_ref(2),
  range_f_ref(2), range_F_ref(2),
  range_f_new(2), range_F_new(2);
  ak.reserve(maxit_ak);
  double err_max;

  range_p[0] = (1.0 - rho) / 2.0;
  range_p[1] = rho + (1.0 - rho) / 2.0;


  double lambda_min = Rcpp::min(lambdas);
  double eps_mell = eps * lambda_min * 0.01;
  double eps_trunc = eps * lambda_min;
  double eps_discr = eps * 0.1;

  NumericVector lambdas_scal = lambdas / lambda_min;

  // compute ak
  if(sum(etas)!=0){
    ak=compute_ak_nc(lambdas_scal,
                         etas,
                         maxit_ak,
                         eps_mell,
                         1.0);
  }else{
    ak=compute_ak_c(lambdas_scal,
                        maxit_ak,
                        eps_mell,
                        1.0);
  }
  if(ak.size() == maxit_ak){// stop with error if the ak did not converged
    stop("Computation of the a_k coefficients did not converged: consider to increase 'maxit_comp'");
  }

  // initial Mellin_ref
  List Mellin_init = Mellin_QF(
    lambdas_scal,  ak,
    maxit_ak,
    h, delta,
    eps_trunc, 1.0);

  // compute function values
  range_q_ref = qQF_c(range_p, Mellin_init, eps_quant, maxit_quant, sum(lambdas_scal));

  List Mellin_ref = Mellin_QF_error(
    lambdas_scal,  ak,
    maxit_ak,
    h, delta,
    eps_trunc, 1.0,
    range_q_ref);
  range_q_ref = qQF_c(range_p, Mellin_ref, eps_quant, maxit_quant, sum(lambdas_scal));
  range_f_ref = dQF_c(range_q_ref, Mellin_ref);
  for(int k=0; k<2; k++) {
    range_f_ref[k] *= 1.0 / lambda_min;
  }

  range_F_ref = pQF_c(range_q_ref, Mellin_ref);


  // first iteration
  delta = std::complex<double> (0.0,step_delta) + delta;
  int direction=1;
  List Mellin_new = Mellin_QF_error(
    lambdas_scal,  ak,
    maxit_ak,
    h, delta,
    eps_trunc, 1.0,
    range_q_ref);
  // new ranges f and F
  range_f_new = dQF_c(range_q_ref, Mellin_new);
  for(int k=0; k<2; k++) {
    range_f_new[k] *= 1.0 / lambda_min;
  }
  range_F_new = pQF_c(range_q_ref, Mellin_new);

  //find maximum error
  err_max = find_maximum_error(range_f_ref, range_F_ref,
                                   range_f_new,range_F_new);

  // decide direction: 0 if the initial value is too high
  if (err_max > eps_discr){
    direction = 0;
  }

  if(direction == 0){// reduce delta
    delta = delta - std::complex<double>(0.0,step_delta);
    for(int i=0; i<maxit_delta; i++){
      delta = delta / 2.0;
      // Compute new quantities
      List Mellin_new = Mellin_QF_error(
        lambdas_scal,  ak,
        maxit_ak,
        h, delta,
        eps_trunc, 1.0,
        range_q_ref);
      range_q_ref = qQF_c(range_p, Mellin_new, eps_quant,
                              maxit_quant,sum(lambdas_scal));
      range_f_new = dQF_c(range_q_ref, Mellin_new);
      for(int k=0; k<2; k++) {
        range_f_new[k] *= 1.0 / lambda_min;
      }

      range_F_new = pQF_c(range_q_ref, Mellin_new);
      range_f_ref = dQF_c(range_q_ref, Mellin_ref);
      for(int k=0; k<2; k++) {
        range_f_ref[k] *= 1.0 / lambda_min;
      }
      range_F_ref = pQF_c(range_q_ref, Mellin_ref);


      //find maximum error
      err_max = find_maximum_error(range_f_ref, range_F_ref,
                                       range_f_new,range_F_new);
      if (err_max < eps_discr){
        break; //stop if the error is lower that the threshold
      }
      if(i == maxit_delta - 1){// stop with error if maxit_delta is reached without the desired error
        stop("Computation of the integration step 'delta' did not converged: consider to increase 'maxit_comp'");
      }
      // new become reference
      range_F_ref = range_F_new;/// DA TOGLIERE
      range_f_ref = range_f_new;
      for(int k=0; k<2; k++) {
        range_f_ref[k] *= 1.0 / lambda_min;
      }

      Mellin_ref = Mellin_new;
    }
    // definitive quantities
    range_q_ref = qQF_c(range_p, Mellin_new, eps_quant,maxit_quant,sum(lambdas_scal));
    delta = delta * 2.0;
  }else{// increase delta
    for(int i=0; i<maxit_delta; i++){
      delta = std::complex<double>(0.0,step_delta) + delta;
      // Compute new quantities
      Mellin_new = Mellin_QF_error(
        lambdas_scal,  ak,
        maxit_ak,
        h, delta,
        eps_trunc, 1.0,
        range_q_ref);
      range_f_new = dQF_c(range_q_ref, Mellin_new);
      for(int k=0; k<2; k++) {
        range_f_new[k] *= 1.0 / lambda_min;
      }

      range_F_new = pQF_c(range_q_ref, Mellin_new);

      //find maximum error
      err_max = find_maximum_error(range_f_ref, range_F_ref,
                                       range_f_new,range_F_new);


      if (err_max > eps_discr){
        break;
      }
      Mellin_ref=Mellin_new;
    }

    delta = delta - std::complex<double>(0.0,step_delta);
  }

  for(int k=0; k<2; k++){
    range_q_ref[k] *= lambda_min;
    }
  Mellin_ref["lambda_min"] = lambda_min;
  return List::create(
    Named("range_q") = range_q_ref,
    Named("Mellin") = Mellin_ref,
    Named("rho") = rho
  );
}

