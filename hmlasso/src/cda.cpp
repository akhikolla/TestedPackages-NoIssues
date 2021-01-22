//' @useDynLib hmlasso
//' @importFrom Rcpp sourceCpp evalCpp

#include <Rcpp.h>
#include <Rinternals.h>
#include <RcppCommon.h>
#include <Rcpp/Benchmark/Timer.h>
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]
#include <boost/numeric/ublas/matrix.hpp>
using namespace Rcpp;

//' calculate covariance matrix
//'
//' @param X design matrix
//' @return covariance matrix
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix covC(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix out(p, p);
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      double total = 0;
      int num = 0;
      for (int k = 0; k < n; k++) {
        if (!R_IsNA(X(k, i)) & !R_IsNA(X(k, j))) {
          total += X(k, i) * X(k, j);
          num += 1;
        }
      }
      out(i, j) = total / num;
    }
  }
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < i; j++) {
      out(i, j) = out(j, i);
    }
  }
  return out;
}

//' soft thresholding function
//'
//' @param z z
//' @param g gamma
//' @return value
//' @keywords internal
//'
// [[Rcpp::export]]
double softThresholdC(double z, double g) {
  if (g < fabs(z)) {
    if(z > 0) {
      return z - g;
    } else {
      return z + g;
    }
  } else {
    return 0;
  }
}

//' update rule function
//'
//' @param u u
//' @param l1 l1
//' @param l2 l2
//' @param v v
//' @return value
//' @keywords internal
//'
// [[Rcpp::export]]
double updateLassoC(double u, double l1, double l2, double v) {
  if (fabs(u) <= l1) {
    return 0;
  } else {
    int s = 0;
    if (u > 0) {
      s = 1;
    } else {
      s = -1;
    }
    return s * (fabs(u) - l1) / (v * (1 + l2));
  }
}

//' Optimize a linear regression model by coordinate descent algorithm using a covariance matrix
//'
//' @param Gamma covariance matrix of explanatory variables
//' @param gamma covariance vector of explanatory and objective variables
//' @param lambda lambda sequence
//' @param warm warm start direction: "lambda" (default) or "delta"
//' @param delta ratio of regularization between l1 and exclusive penalty terms
//' @param R matrix using exclusive penalty term
//' @param maxit max iteration
//' @param eps convergence threshold for optimization
//' @param init_beta initial values of beta
//' @param strong whether use strong screening or not
//' @return standardized beta
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericMatrix covCdaC(NumericMatrix Gamma, NumericVector gamma, NumericVector lambda,
                      NumericMatrix R, NumericMatrix init_beta, double delta = 0,
                      double maxit = 1e+4, double eps = 1e-04,
                      String warm = "lambda", bool strong = true) {
  // initialize
  int p = Gamma.ncol();
  int nlambda = lambda.size();
  double z, g, resid;
  NumericMatrix beta(p, nlambda);
  NumericVector beta_prev(p);
  LogicalVector candid(p), candid2(p);

  // return beta;

  // CDA for all lambda
  for (int k = 0; k < nlambda; k++) {
    double prod;
    NumericVector beta_now(p);

    // Rcout << "k: " << k << std::endl;

    // beta-init or warm-start
    if (warm=="delta") {
      beta(_, k) = init_beta(_, k);
    } else if (warm=="lambda") {
      if (k != 0) {
        beta(_, k) = beta(_, k - 1);
      }
    }

    // strong screening
    if (strong && (k != 0)) {
      for (int j = 0; j < p; j++) {
        prod = 0;
        for (int l = 0; l < p; l++) {
          prod += Gamma(j, l) * beta(j, k-1);
        }
        candid(j) = (fabs(gamma(j) - prod) >= 2 * lambda(k) - lambda(k-1));
      }
    } else {
      for (int j = 0; j < p; j++) {
        candid(j) = true;
      }
    }

    // CDA
    for (int i = 0; i < maxit; i++) {
      for (int j = 0; j < p; j++) {
        if (beta(j, k)==0) {
          continue;
        }
        prod = 0;
        for (int l = 0; l < p; l++) {
          if ((beta(l, k) != 0) && (l != j)) {
            prod += Gamma(j, l) * beta(l, k);
          }
        }
        z = gamma(j) - prod;
        if (delta == 0) {
          g = lambda(k);
          beta(j, k) = softThresholdC(z, g) / Gamma(j, j);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * fabs(beta(l, k));
            }
          }
          g = (1 + delta * prod) * lambda(k);
          beta(j, k) = softThresholdC(z, g) / (Gamma(j, j) + delta * lambda(k) * R(j, j));
        }
      }
      resid = max(abs(beta(_, k) - beta_prev));
      beta_prev = beta(_, k);
      if (resid < eps) {
        // Rcout << "1 iteration " << i << std::endl;
        break;
      }
    }
    // CDA
    for (int i = 0; i < maxit; i++) {
      for (int j = 0; j < p; j++) {
        if (!candid(j)) {
          continue;
        }
        prod = 0;
        for (int l = 0; l < p; l++) {
          if ((beta(l, k) != 0) && (l != j)) {
            prod += Gamma(j, l) * beta(l, k);
          }
        }
        z = gamma(j) - prod;
        if (delta == 0) {
          g = lambda(k);
          beta(j, k) = softThresholdC(z, g) / Gamma(j, j);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * fabs(beta(l, k));
            }
          }
          g = (1 + delta * prod) * lambda(k);
          beta(j, k) = softThresholdC(z, g) / (Gamma(j, j) + delta * lambda(k) * R(j, j));
        }
      }
      resid = max(abs(beta(_, k) - beta_prev));
      beta_prev = beta(_, k);
      if (resid < eps) {
        // Rcout << "2 iteration " << i << std::endl;
        break;
      }
    }
    // CDA
    for (int i = 0; i < maxit; i++) {
      for (int j = 0; j < p; j++) {
        prod = 0;
        for (int l = 0; l < p; l++) {
          if ((beta(l, k) != 0) && (l != j)) {
            prod += Gamma(j, l) * beta(l, k);
          }
        }
        z = gamma(j) - prod;
        if (delta == 0) {
          g = lambda(k);
          beta(j, k) = softThresholdC(z, g) / Gamma(j, j);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * fabs(beta(l, k));
            }
          }
          g = (1 + delta * prod) * lambda(k);
          beta(j, k) = softThresholdC(z, g) / (Gamma(j, j) + delta * lambda(k) * R(j, j));
        }
      }
      resid = max(abs(beta(_, k) - beta_prev));
      beta_prev = beta(_, k);
      if (resid < eps) {
        // Rcout << "3 iteration " << i << std::endl;
        break;
      }
    }

  }

  return beta;
}
