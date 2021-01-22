#include <RcppArmadillo.h>

// [[Rcpp::export(name = ".estimate_VECM")]]
SEXP estimate_VECM(arma::mat y, 
                   int r,
                   int p, // the number of lags
                   int model, 
                   arma::mat dummy) {
  
    int n = y.n_rows, k = y.n_cols, T = n - p, n_dummy = dummy.n_cols;
  
  // sets number of collumns for z1 and z2:
  int z1_cols, z2_cols; 
  if(model == 1) z1_cols = k; 
  else z1_cols = k + 1;
  if(model == 3) z2_cols = k * (p - 1) + 1 + n_dummy; 
  else z2_cols = k * (p - 1) + n_dummy;
  
  arma::mat z0(T, k), z1(T, z1_cols), z2(T, z2_cols);
  
  // computes the difference matrix of y:
  arma::mat y_diff = y.submat(1, 0, n - 1, k-1) - y.submat(0, 0, n - 2, k - 1);
  
  // populates the z-matrices:
  z0 = y_diff(arma::span(p - 1, n - 2), arma::span(0, k - 1));
  if(model == 1){
    z1(arma::span(0, T - 1), arma::span(0, k - 1)) = y(arma::span(p - 1, n - 2), arma::span(0, k - 1));

    for(int j = 0; j < p - 1; ++j){
      z2(arma::span(0, T - 1), arma::span(j * k, (j + 1) * k - 1)) =
        y_diff(arma::span(p - j - 2, n - 3 - j), arma::span(0, k - 1));
    }
  } else if(model == 2){
    z1(arma::span(0, T - 1), arma::span(0, k - 1)) = y(arma::span(p - 1, n - 2), arma::span(0, k - 1));
    z1.col(k) = arma::ones(T); //constant
    
    for(int j = 0; j < p - 1; ++j){
      z2(arma::span(0, T - 1), arma::span(j * k, (j + 1) * k - 1)) =
        y_diff(arma::span(p - j - 2, n - 3 - j), arma::span(0, k - 1));
    }
  } else if(model == 3){
    z1(arma::span(0, T - 1), arma::span(0, k - 1)) = y(arma::span(p - 1, n - 2), arma::span(0, k - 1));
    z1.col(k) = arma::regspace(p + 1, n); // trend
    
    for(int j = 0; j < p - 1; ++j){
      z2(arma::span(0, T - 1), arma::span(j * k, (j + 1) * k - 1)) =
        y_diff(arma::span(p - j - 2, n - 3 - j), arma::span(0, k - 1));
    }
    z2.col((p - 1) * k) = arma::ones(T);  //constant
  }
  
  // possibly adds dummy variables:
  if(n_dummy > 0) z2.cols(z2_cols - n_dummy, z2_cols - 1) = 
    dummy(arma::span(p, n - 1), arma::span(0, n_dummy - 1));
  
  // computes the m**, r* and s** matrices:
  arma::mat m02 = z0.t() * z2 / T;
  arma::mat m12 = z1.t() * z2 / T;
  arma::mat m22 = z2.t() * z2 / T;
  arma::mat r0 = z0 - z2 * (m02 * arma::inv_sympd(m22)).t();
  arma::mat r1 = z1 - z2 * (m12 * arma::inv_sympd(m22)).t();
  arma::mat s00 = r0.t() * r0 / T;
  arma::mat s01 = r0.t() * r1 / T;
  arma::mat s10 = s01.t();
  arma::mat s11 = r1.t() * r1 / T;
  
  // computes the generalized eigen values and vectors:
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_pair(eigval, eigvec, s10 * arma::inv_sympd(s00) * s01, s11);
  
  // sorts the eigen values and vectors, and transforms to non-complex values:
  arma::uvec indices = sort_index(eigval, "descend");
  arma::vec eigval_dbl = arma::conv_to<arma::mat>::from(eigval(indices));
  arma::mat eigvec_dbl = arma::conv_to<arma::mat>::from(eigvec.cols(indices));
  
  // normalizes the eigenvector by the diagonal:
  eigvec_dbl = eigvec_dbl * arma::inv(arma::diagmat(eigvec_dbl));
   
  // extracts 'beta' (including 'rho' part) from the eigen vectors:
  arma::mat beta_and_rho;
  if(r > 0){
    beta_and_rho = eigvec_dbl.cols(arma::span(0, r - 1));
  }
  
  // computes 'alpha' and 'psi':
  arma::mat alpha;
  arma::mat psi;
  if(r > 0){
    alpha = s01 * beta_and_rho * inv(beta_and_rho.t() * s11 * beta_and_rho);
    psi = m02 * inv(m22) - alpha * beta_and_rho.t() * m12 * inv(m22);
  } else{
    psi = m02 * inv(m22);
  }
  
  // splits 'beta_and_rho' into 'beta' and 'rho':
  arma::mat beta;
  arma::mat rho;
  if(r > 0){
    if(model == 2 || model == 3){
      rho = beta_and_rho.row(beta_and_rho.n_rows - 1);
      beta = beta_and_rho.rows(arma::span(0, beta_and_rho.n_rows - 2));
    } else if(model == 1){
      beta = beta_and_rho.rows(arma::span(0, beta_and_rho.n_rows - 1));
    }
  }
  
  // splits 'psi' into 'gamma' and 'phi':
  arma::mat gamma;
  arma::mat phi;
  if(p > 1){
    gamma = psi.cols(arma::span(0, (p - 1) * k - 1));
    if(model == 3){
      phi = psi.col((p - 1) * k);
    } 
  }
  
  arma::mat dummy_phi; //the matrix of dummy parameters
  if(n_dummy > 0){
    //dummy_phi = psi.cols(arma::span((p - 1) * k + 1, psi.n_cols - 1));
    dummy_phi = psi.cols(arma::span(psi.n_cols - n_dummy, psi.n_cols - 1));
  } 
  
  // computes the residauls:
  arma::mat residuals;
  if(r > 0){
    residuals = z0 - z1 * beta_and_rho * alpha.t()  - z2 * psi.t();
  } else{
    residuals = z0 - z2 * psi.t();
  }
  
  // conputes the Q statistic:
  double Q = 0;
  if(r < k){
    Q = - T * arma::accu(log(arma::ones(k - r) - eigval_dbl(arma::span(r, k - 1))));
  }
  
  Rcpp::List returnList;

  returnList["eigval_dbl"] =  eigval_dbl;
  returnList["eigvec_dbl"] =  eigvec_dbl;
  returnList["beta"] =  beta;
  returnList["rho"] =  rho;
  returnList["alpha"] =  alpha;
  returnList["Q"] =  Q;
  returnList["gamma"] =  gamma;
  returnList["phi"] =  phi;
  returnList["dummy_coefs"] =  dummy_phi;
  returnList["residuals"] =  residuals;
  returnList["psi"] =  psi;
  
  return(returnList);
}

// [[Rcpp::export(name = ".make_VECM")]]
arma::mat make_VECM(arma::mat Ystart,
               arma::mat e,
               arma::mat alpha,
               arma::mat beta,
               arma::mat gamma,
               arma::mat rho,
               arma::mat phi,
               int model){


  int k = e.n_cols, p = gamma.n_cols / k + 1, n = e.n_rows + p, r = beta.n_cols;

  arma::mat y(n, k);
  arma::mat delta_y(n, k);

  // fills the first rows of 'y' with 'Ystart':
  y(arma::span(0, p - 1), arma::span(0, k - 1)) = Ystart;
  
  // fills the first rows of 'delta_y' with delta 'Ystart':
  delta_y.rows(arma::span(1, p - 1)) =
    Ystart.rows(arma::span(1, p - 1)) - Ystart.rows(arma::span(0, p - 2));
  
  // adds the errors:
  delta_y(arma::span(p, n - 1), arma::span(0, k - 1)) = e;
  
  // adds constants:
  if(model == 2 && r > 0){
    arma::mat alpha_rho_t = rho * alpha.t();
    for(int i = p - 1; i < n; ++i) delta_y.row(i) = delta_y.row(i) + alpha_rho_t;
  }
  if(model == 3){
    for(int i = p - 1; i < n; ++i) delta_y.row(i) = delta_y.row(i) + phi.t();
  }
  
  // adds trend:
  if(model == 3 && r > 0){
    arma::mat alpha_rho_t = rho * alpha.t();
    for(int i = p - 1; i < n; ++i) y.row(i) = y.row(i) + alpha_rho_t * i;
  }
  
  // adds the dynamic part:
  arma::mat pi_t;
  if(r > 0) pi_t = beta * alpha.t();
  for(int i = p; i < n; ++i){
    
    if(r > 0) delta_y.row(i) += y.row(i - 1) * pi_t; 
    
    for(int j = 0; j < p - 1; ++j){
      delta_y.row(i) += delta_y.row(i - j - 1) * gamma.cols(arma::span(j * k, (j + 1) * k - 1)).t(); 
    }
    
    y.row(i) = y.row(i - 1) + delta_y.row(i);
  }
  
  return(y);
}

// [[Rcpp::export(name = ".vecm_check_eigen")]]
SEXP vecm_check_eigen(arma::mat alpha,
                      arma::mat beta,
                      arma::mat gamma,
                      int r,
                      int k,
                      int p){
  
  arma::mat A(p * k, p * k, arma::fill::zeros);
  
  A(arma::span(0, k - 1), arma::span(0, k - 1)) = arma::eye(k, k) + 
    gamma(arma::span(0, k - 1), arma::span(0, k - 1));
  
  if(r > 0) A(arma::span(0, k - 1), arma::span(0, k - 1)) = A(arma::span(0, k - 1), arma::span(0, k - 1)) + 
    alpha * beta.t();
  
  if(p > 2){ 
    for(int j = 1; j < p - 1; ++j){
      A(arma::span(0, k - 1), arma::span(k * j, k * (j + 1) - 1)) = 
        gamma(arma::span(0, k - 1), arma::span(k * j, k * (j + 1) - 1)) -
        gamma(arma::span(0, k - 1), arma::span(k * (j - 1), k * j - 1));
    }
  }
  
  if(p > 1){
    A(arma::span(0, k - 1), arma::span(k * (p - 1), k * p - 1)) =
      -gamma(arma::span(0, k - 1), arma::span(k * (p - 2), k * (p - 1) - 1));
      
    A(arma::span(k, p * k - 1), arma::span(0, k * (p - 1) - 1)) = arma::eye(k * (p - 1), k * (p - 1));
  }
  
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  
  arma::eig_gen(eigval, eigvec, A);
  
  Rcpp::List returnList;
  
  returnList["eigval"] =  eigval;
  returnList["abs_eigval"] =  abs(eigval);
  returnList["eigvec"] =  eigvec;
  returnList["A"] =  A;
  
  return(returnList);
}
