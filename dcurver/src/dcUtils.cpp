# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

// using namespace Rcpp;

arma::mat expVec(double x, int deg) {
  // Create a column vector of the exponential series:
  // 1 + theta + theta^2 + ... + theta_k.
  
  arma::mat out(deg+1, 1);
  
  for (int i = 0; i < deg+1; i++) {
    out(i, 0) = pow(x, i);
  }
  
  return out;
}

arma::mat cMat (int k, NumericVector phi) {
  // Create c matrix.
  
  // c is a term in Equation (9) of Woods & Lin (2009).
  // This function creates the c matrix given a degree of k and a vector phi of
  // length k, containing the parameters of the polynomial P_k^2(\theta).
  arma::mat prim_c(1,1);
  arma::mat c(1, 1);
  
  prim_c[1] = 1;
  
  // Fill primitive prim_c matrix
  if (k != 0) {
    prim_c = arma::mat(k+1, k);
    c = arma::mat(k+1, 1);
    
    for (unsigned int i = 0; i < prim_c.n_rows; i++) {
      if (i + 1 == prim_c.n_rows) { // if last row
        // write all cos's
        for (unsigned int j = 0; j < prim_c.n_cols; j++) {
          prim_c(i, j) = cos(phi[j]);
        }
      } else { // if not last row
        for (unsigned int j = 0; j < prim_c.n_cols; j++) {
          if (j > i) {
            prim_c(i, j) = 1;
          } else if (j == i) { // if last col
            prim_c(i, j) = sin(phi[j]);
          } else { // if not last col
            prim_c(i, j) = cos(phi[j]);
          }
        }
      }
    }
  }
  
  // Collapse cols by multiplication to create c matrix
  for (unsigned int i = 0; i < prim_c.n_rows; i++) {
    double mult = 1;
    
    for (unsigned int j = 0; j < prim_c.n_cols; j++) {
      mult = mult * prim_c(i, j);
    }
    
    c(i, 0) = mult;
  }
  
  return c;
}

arma::vec insZ (arma::vec vals) {
  arma::vec tmp(vals.n_elem * 2, arma::fill::zeros);
  
  for (unsigned int i = 0; i < vals.n_elem; i++) {
    tmp[i*2] = vals(i);
  }
  
  return tmp;
}

arma::mat create_M (int k) {
  // Function to create M given the number of parameters.
  
  // M is defined in Equation 9 in Woods & Lin (2009).
  
  arma::mat M(k+1, k+1, arma::fill::zeros);
  arma::rowvec valsZ;
  arma::vec vals(k+1);

  // M matrix values.
  if (k == 1) {
    vals = arma::vec({1, 1}); // E[X^0, X^2]
  } else if (k == 2) {
    vals = arma::vec({1, 1, 3}); // E[X^0], E[X^2], E[X^4]
  } else if (k == 3) {
    vals = arma::vec({1, 1, 3, 15}); // ... and so on.
  } else if (k == 4) {
    vals = arma::vec({1, 1, 3, 15, 105});
  } else if (k == 5) {
    vals = arma::vec({1, 1, 3, 15, 105, 945});
  } else if (k == 6) {
    vals = arma::vec({1, 1, 3, 15, 105, 945, 10395});
  } else if (k == 7) {
    vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135});
  } else if (k == 8) {
    vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135, 2027025});
  } else if (k == 9) {
    vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425});
  } else if (k == 10) {
    vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425, 654729075});
  }
  
  valsZ = insZ(vals).t();
  
  for (int i = 0; i < k+1; i++) {
    M.row(i) = valsZ.head(k+1);
    valsZ = arma::shift(valsZ, -1);
  }
  
  return M;
}

arma::mat invBMat (int k) {
  // Obtain B matrix and find the inverse.
  
  // Following the notation in Woods & Lin (2009).
  // Inverse of B is a frequently used term. B itself is defined,
  // B^T B = M. M is given by create_M() defined above. Then, B is obtained
  // by Cholesky decomposition and inversed.
  arma::mat B(k+1, k+1);
  arma::mat invB(k+1, k+1);
  arma::mat M(k+1, k+1);
  
  M = create_M(k);
  B = arma::chol(M);
  invB = arma::inv(B);
  
  return invB;
}

