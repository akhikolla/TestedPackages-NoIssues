// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matOmega1(arma::mat W, double rho) {
  // Omega1 - SAR(1)
  int n = W.n_rows;
  arma::mat I = arma::eye<arma::mat>(n, n);
  arma::mat S1 = I - rho * W;
  return arma::inv_sympd(S1.t() * S1);
}

// [[Rcpp::export]]
arma::mat matOmega2(int nTime, double rho) {
  // Omega2 - AR(1)
  arma::mat Ome2(nTime, nTime);
  Ome2.fill(0.0);
  for(int i = 0; i < nTime; ++i) {
    Ome2.diag(i) += pow(rho, i);
  }
  Ome2 += Ome2.t();
  Ome2.diag() *= 0.0;
  Ome2.diag() += 1;
  return 1/(1-pow(rho, 2)) * Ome2;
}

// [[Rcpp::export]]
arma::mat matBlockDiagonal(arma::mat X, int n) {
  arma::mat XX(X.n_cols*n, X.n_cols*n);
  XX.fill(0.0);
  for(int r = 0; r < n; ++r) {
    XX.submat(arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1), arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1)) =
      X;
  }
  return XX;
}

// [[Rcpp::export]]
Rcpp::List matVInvT(arma::mat Ome1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  // Inverse of variance covariance matrix of ?-temporal model
  // This can work for the temporal model if Ome1 is an identity and the spatial
  // model when Ome1 is the coresponding variance matrix
  arma::mat Ome2 = matOmega2(Z1.n_rows / Ome1.n_rows, rho2);
  arma::mat Ad = sigma2 * Ome2;
  arma::mat A = matBlockDiagonal(Ad, Ome1.n_rows);
  A.diag() += sigmaSamplingError;
  arma::mat Ainv = A;

  for(int r = 0; r < Ome1.n_rows; ++r) {
    Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)) =
      arma::inv(Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)));
  }

  arma::mat V = sigma1 * Z1 * Ome1 * Z1.t() + A;
  arma::mat AinvZ1 = Ainv * Z1;
  arma::mat Ome1inv = arma::inv(sigma1 * Ome1);
  arma::mat Vinv = Ainv - AinvZ1 * arma::inv(Ome1inv + Z1.t() * AinvZ1) * AinvZ1.t();
  return Rcpp::List::create(Rcpp::Named("V", V),
                            Rcpp::Named("VInv", Vinv));
}

// [[Rcpp::export]]
arma::mat matVDerS1(arma::mat Ome1, arma::mat Z1) {
  // derivative of V with respect to sigma1
  return Z1 * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVDerS2(arma::mat Ome2, int nDomains) {
  // derivative of V with respect to sigma2
  return matBlockDiagonal(Ome2, nDomains);
}

// [[Rcpp::export]]
arma::mat matVDerR1(double rho1, double sigma1, arma::mat Z1, arma::mat Ome1, arma::mat W) {
  // derivative of V with respect to rho1
  return -sigma1 * Z1 * Ome1 * (-W-W.t() + 2 * rho1 * W.t() * W) * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVDerR2(double rho2, double sigma2, arma::mat Ome2, int nDomains) {
  // derivative of V with respect to rho2
  arma::mat ome2derR2(Ome2.n_cols, Ome2.n_cols);
  ome2derR2.fill(0.0);
  for(int i = 1; i < Ome2.n_cols; ++i) {
    ome2derR2.diag(i) += i * pow(rho2, i-1);
  }

  ome2derR2 += ome2derR2.t();
  ome2derR2 = 1/(1-pow(rho2, 2)) * (ome2derR2 + 2 * rho2 * Ome2);
  return sigma2 * matBlockDiagonal(ome2derR2, nDomains);
}

arma::colvec psiOne(arma::colvec u, double k = 1.345){
  //arma::colvec weights = rep(1, length(u));
  //sm<-median(abs(u/sqrt(var.weights)))/0.6745
  Rcpp::NumericVector tmp = as<NumericVector>(wrap(k/abs(u)));
  arma::colvec w = Rcpp::pmin(1.0, tmp);
  return w%u;
}

// [[Rcpp::export]]
arma::colvec fixedPointSigma(
    arma::colvec y, arma::mat X, arma::colvec beta, arma::colvec sigma, double rho,
    arma::mat Z1, arma::mat Ome1, arma::mat Z, arma::colvec sigmaSamplingError,
    double k, double K) {

  // Variance-Covarianze matrix
  Rcpp::List Vlist = matVInvT(Ome1, sigma(0), rho, sigma(1), Z1, sigmaSamplingError);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);

  // sqrt of U + inverse
  arma::mat sqrtU(V.n_rows, V.n_rows);
  sqrtU.fill(0.0);
  sqrtU.diag() = sqrt(V.diag());

  arma::mat sqrtUinv(V.n_rows, V.n_rows);
  sqrtUinv.fill(0.0);
  sqrtUinv.diag() = 1/sqrtU.diag();

  // residuals and huber
  arma::colvec resid = sqrtUinv * (y - X * beta);
  arma::colvec phiR = psiOne(resid);

  // Derivatives
  arma::mat Ome2 = matOmega2(Z1.n_rows / Ome1.n_rows, rho);
  arma::mat derVSigma1 = matVDerS1(Ome1, Z1);
  arma::mat derVSigma2 = matVDerS2(Ome2, Ome1.n_rows);

  // ZVuZt + inverse
  arma::mat ZVuZt = V;
  ZVuZt.diag() -= sigmaSamplingError;
  arma::mat ZVuZtinv = arma::inv(trimatu(chol(ZVuZt)));
  ZVuZtinv = ZVuZtinv * ZVuZtinv.t();

  // OmegaBar
  arma::mat Ome1Bar(Z.n_cols, Z.n_cols);
  arma::mat Ome2Bar(Z.n_cols, Z.n_cols);
  Ome1Bar.fill(0.0);
  Ome2Bar.fill(0.0);

  Ome1Bar(arma::span(0, Ome1.n_rows-1), arma::span(0, Ome1.n_rows-1)) = Ome1;
  Ome2Bar(arma::span(Ome1.n_rows, Ome2Bar.n_rows-1), arma::span(Ome1.n_rows, Ome2Bar.n_rows-1)) = matBlockDiagonal(Ome2, Ome1.n_rows);

  // Compute a(theta) and A(theta)
  arma::mat tmp1 = phiR.t() * sqrtU * Vinv;

  arma::colvec a(2);

  a(0) = as_scalar(tmp1 * derVSigma1 * tmp1.t());
  a(1) = as_scalar(tmp1 * derVSigma2 * tmp1.t());

  arma::mat tmp2 = K * Vinv * derVSigma1 * ZVuZtinv;
  arma::mat tmp3 = K * Vinv * derVSigma2 * ZVuZtinv;

  arma::mat tmp4 = Z * Ome1Bar * Z.t();
  arma::mat tmp5 = Z * Ome2Bar * Z.t();

  arma::mat A(2,2);
  A(0,0) = trace(tmp2 * tmp4);
  A(0,1) = trace(tmp2 * tmp5);
  A(1,0) = trace(tmp3 * tmp4);
  A(1,1) = trace(tmp3 * tmp5);

  return inv(A) * a;

}
