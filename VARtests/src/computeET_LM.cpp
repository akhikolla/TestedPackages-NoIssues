#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat computeET_LM(arma::mat resid, int h) {
  
  int n = resid.n_rows, K = resid.n_cols, p = (h + 1) * K + K * (K - 1) / 2;
  
  arma::mat score(p, 1, arma::fill::zeros);
  arma::mat info(p, p, arma::fill::zeros);
  arma::mat dVecDinvdTheta(K * K, p, arma::fill::zeros);
  arma::mat dvecPdTheta(K * K, p, arma::fill::zeros);
  
  arma::mat cov = arma::cov(resid);
  
  arma::mat resid2 = arma::fliplr(pow(resid, 2).t());
  arma::mat ee(K, K);
  
  arma::mat D = sqrt(diagmat(cov));
  arma::mat P = arma::cor(resid);
  
  arma::mat Pinv = arma::inv_sympd(P);
  arma::mat Dinv = arma::inv(diagmat(D));
  
  arma::mat DinvPinv = Pinv * Dinv;
  arma::mat PinvDinv = Dinv * Pinv;
  
  int i = 0, j = 0, row = 0, colTemp = 0;
  
  for(j = 0; j < K; ++j) {
    for(i = j + 1; i < K; ++i) {
      colTemp = K * j - (j * j + 3 * j) / 2 + i + (h + 1) * K - 1;
      dvecPdTheta(j * K + i, colTemp) = 1;
      dvecPdTheta(i * K + j, colTemp) = 1;
    }
  }
  
  arma::mat factor1 = arma::kron(D, D) + .5 *(arma::kron(Pinv, cov) + arma::kron(cov, Pinv));
  arma::mat factor2 = .5 * (arma::kron(D, Pinv) + arma::kron(Pinv, D));
  arma::mat factor3 = .5 * dvecPdTheta.t() * arma::kron(Pinv, Pinv) * dvecPdTheta;
  
  for(row = h; row < n; ++row) {
    
    // builds the dVecDinvdTheta matrix:
    for(i = 0; i < K; ++i) {
      
      dVecDinvdTheta(i*(K + 1), i * (h + 1)) = - 1 / (2 * pow(D(i, i), 3));
      dVecDinvdTheta(i*(K + 1), arma::span(i * (h + 1) + 1, (i + 1)*(h + 1) - 1)) =
        - resid2(i, arma::span(n - row, n - row + h - 1)) / (2 * pow(D(i, i), 3));
    }
    
    // outer product of the residuals for row 'row':
    ee = resid.row(row).t() * resid.row(row);    
    
    // score vector:
    score = score 
            + dVecDinvdTheta.t() * vectorise(D - .5 * ee * DinvPinv - .5 * PinvDinv * ee)
            + .5 * dvecPdTheta.t() * vectorise(PinvDinv * ee * DinvPinv - Pinv);
    
    // information matrix:
    info = info
          + dVecDinvdTheta.t() * factor1 * dVecDinvdTheta
          - dVecDinvdTheta.t() * factor2 * dvecPdTheta
          - dvecPdTheta.t() * factor2 * dVecDinvdTheta
          + factor3;
  }
  
  score = score / n;
  info = info / n;
  
  return n * score.t() * arma::inv_sympd(info) * score;
}
