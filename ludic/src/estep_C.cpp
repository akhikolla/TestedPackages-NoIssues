#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//'Fast C++ implementation of Winkler's Method E step
//'
//'@keywords internal
// [[Rcpp::export]]
arma::mat estep_C_vect(arma::mat agreemat, double p, arma::colvec m, arma::colvec u){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  mat res(N, 2);
  
  mat agreemat_neg = 1.0 - agreemat;
  vec a = exp(log(p) + agreemat*log(m) + agreemat_neg*log(1.0-m));
  vec b = exp(log(1.0-p) + agreemat*log(u) + agreemat_neg*log(1.0-u));
  res.col(0) = a/(a+b);
  res.col(1) = b/(a+b);
  
  return(res);
}



//'C++ implementation of the E and M steps from Winkler's EM algorithm estimating FS method 
//'using sparse matrices for big sample sizes
//'
//'@keywords internal
// [[Rcpp::export]]

List EMstep_C_sparse_big(arma::mat mat_A, arma::mat mat_B, double p, arma::rowvec m, arma::rowvec u){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  vec g_m(N);
  vec g_u(N);
  
  //E step
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      rowvec agreerow_pairl_neg(abs(mat_A.row(i) - mat_B.row(j)));
      rowvec agreerow_pairl(1.0 - agreerow_pairl_neg);
      double a = exp(log(p) + sum(agreerow_pairl%log(m) + (agreerow_pairl_neg)%log(1.0 - m)));
      double b = exp(log(1.0 - p) + sum(agreerow_pairl%log(u) + (agreerow_pairl_neg)%log(1.0 - u)));
      g_m(l) = a/(a+b);
      g_u(l) = b/(a+b);
    }
  }
  
  vec m_res(zeros(K));
  vec u_res(zeros(K));
  
  //M step
  double gm_summed = sum(g_m);
  double gu_summed = sum(g_u);
  double p_res = gm_summed/N;
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      colvec agreerow_pairl(conv_to<colvec>::from(1.0 - abs(mat_A.row(i) - mat_B.row(j))));
      m_res += g_m(l)*agreerow_pairl / gm_summed;
      u_res += g_u(l)*agreerow_pairl / gu_summed;
    }
  }
  return(Rcpp::List::create(Rcpp::Named("p") = p_res, Rcpp::Named("m") = m_res.t(), 
                            Rcpp::Named("u") = u_res.t()));
}

