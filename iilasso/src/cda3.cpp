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
//' 
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
          beta(j, k) = softThresholdC(z, g);
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
          beta(j, k) = softThresholdC(z, g);
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
          beta(j, k) = softThresholdC(z, g);
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


//' Optimize a logistic regression model by coordinate descent algorithm using a design matrix
//' 
//' @param X_tilde standardized matrix of explanatory variables
//' @param y vector of objective variable
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
//' 
// [[Rcpp::export]]
NumericMatrix logitCdaC(NumericMatrix X_tilde, NumericVector y, NumericVector lambda,
                        NumericMatrix R, NumericMatrix init_beta, double delta = 0,
                        double maxit = 1e+4, double eps = 1e-04,
                        String warm = "lambda", bool strong = true) {
  // initialize
  int n = X_tilde.nrow();
  int p = X_tilde.ncol() + 1;
  int nlambda = lambda.size();

  // X
  NumericMatrix X_tilde2(n, p);
  for (int i = 0; i < n; i++) {
    X_tilde2(i, 0) = 1;
  }
  for (int j = 1; j < p; j++) {
    for (int i = 0; i < n; i++) {
      X_tilde2(i, j) = X_tilde(i, j - 1);
    }
  }

  // beta
  NumericMatrix beta(p, nlambda);
  NumericVector beta_tilde(p);
  NumericVector eta(n);

  // similarity matrix
  NumericMatrix R2(p, p);
  for (int j = 0; j < p; j++) {
    R2(0, j) = 0;
  }
  for (int j = 0; j < p; j++) {
    R2(0, j) = 0;
  }
  for (int j1 = 1; j1 < p; j1++) {
    for (int j2 = 1; j2 < p; j2++) {
      R2(j1, j2) = R(j1 - 1, j2 - 1);
    }
  }

  // candidate variables
  LogicalVector candid(p);
  candid(0) = true;
  NumericVector pi(n), w(n), s(n), r(n);

  // CDA for all lambda
  for (int k = 0; k < nlambda; k++) {
    // Rcout << "lambda " << k << std::endl;

    // beta-init or warm-start
    if (warm=="delta") {
      beta(_, k) = init_beta(_, k);
    } else if (warm=="lambda") {
      if (k != 0) {
        beta(_, k) = beta(_, k - 1);
      }
    }

    int iter = 0;
    while(iter < maxit) {
      while(iter < maxit) {
        iter++;

        beta_tilde = beta(_, k);

        for (int i = 0; i < n; i++) {
          eta(i) = 0;
          for (int j = 0; j < p; j++) {
            eta(i) = eta(i) + X_tilde2(i, j) * beta_tilde(j);
          }
          if (eta(i) > 10) {
            pi(i) = 1;
          } else if (eta(i) < -10) {
            pi(i) = 0;
          } else {
            pi(i) = 1 / (1 + exp(- eta(i)));
          }
          w(i) = pi(i) * (1 - pi(i));
          if (w(i) < 0.0001) {w(i) = 0.0001;}
          s(i) = y(i) - pi(i);
          r(i) = s(i) / w(i);
        }

        double maxChange = 0;

        for (int j = 0; j < p; j++) {
          if(candid(j)) {
            double xwr = 0;
            for (int i = 0; i < n; i++) {
              xwr = xwr + X_tilde2(i, j) * w(i) * r(i);
            }
            double xwx = 0;
            for (int i = 0; i < n; i++) {
              xwx = xwx + X_tilde2(i, j) * w(i) * X_tilde2(i, j);
            }

            double u, v, l1, l2;
            if (delta == 0) { // normal lasso
              u = xwr / n + (xwx / n) * beta_tilde(j);
              v = xwx / n;
              l1 = lambda(k);
              l2 = 0;
            } else { // exclusive lasso
              u = xwr / n + (xwx / n) * beta_tilde(j);
              v = xwx / n + delta * lambda(k) * R2(j, j);
              double rb = 0;
              for (int jj = 0; jj < p; jj++) {
                if ((jj != j) & (beta(jj, k) != 0)) {
                  rb = rb + R2(j, jj) * fabs(beta(jj, k));
                }
              }
              l1 = lambda(k) + delta * lambda(k) * rb;
              l2 = 0;
            }

            if (j == 0) {
              // beta(j, k) = updateLassoC(u, 0, l2, 1);
              // beta(j, k) = u / v;
              beta(j, k) = xwr / xwx + beta_tilde(j);
            } else {
              beta(j, k) = updateLassoC(u, l1, l2, v);
            }

            double shift = beta(j, k) - beta_tilde(j);
            if (shift != 0) {
              for (int i = 0; i < n; i++) {
                double si = shift * X_tilde2(i, j);
                r(i) = r(i) - si;
                eta(i) = eta(i) + si;
              }
            }

            if (maxChange < fabs(shift) * sqrt(v)) {
              maxChange = fabs(shift) * sqrt(v);
            }
          }
        }

        // check for convergence
        // double maxChange = max(abs(beta(_, k) - beta_tilde));
        beta_tilde = beta(_, k);
        if (maxChange < eps) {
          // Rcout << iter << std::endl;
          break;
        }
      }

      int violations = 0;
      for (int j = 0; j < p; j++) {
        if (!candid(j)) {
          double z = 0;
          for (int i = 0; i < n; i++) {
            z = z + X_tilde2(i, j) * s(i) / n;
          }
          if (fabs(z) > lambda(k)) {
            candid(j) = true;
            violations++;
          }
        }
      }
      if (violations == 0) {
        // Rcout << iter << std::endl;
        break;
      }
    }

    // break if maxit reached
    if (iter >= maxit) {
      return(beta);
    }
  }

  return beta;

}

//' (Experimental) Optimize an ULasso linear regression problem by coordinate descent algorithm using a covariance matrix
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
//' 
// [[Rcpp::export]]
NumericMatrix covCdaC2(NumericMatrix Gamma, NumericVector gamma, NumericVector lambda,
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
          beta(j, k) = softThresholdC(z, g);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * beta(l, k);
            }
          }
          z = z - delta * lambda(k) * prod;
          g = lambda(k);
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
          beta(j, k) = softThresholdC(z, g);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * beta(l, k);
            }
          }
          z = z - delta * lambda(k) * prod;
          g = lambda(k);
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
          beta(j, k) = softThresholdC(z, g);
        } else {
          prod = 0;
          for (int l = 0; l < p; l++) {
            if ((beta(l, k) != 0) && (l != j)) {
              prod += R(j, l) * beta(l, k);
            }
          }
          z = z - delta * lambda(k) * prod;
          g = lambda(k);
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


//' (Experimental) Optimize an ULasso logistic regression problem by coordinate descent algorithm using a design matrix
//' 
//' @param X_tilde standardized matrix of explanatory variables
//' @param y vector of objective variable
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
//' 
// [[Rcpp::export]]
NumericMatrix logitCdaC2(NumericMatrix X_tilde, NumericVector y, NumericVector lambda,
                        NumericMatrix R, NumericMatrix init_beta, double delta = 0,
                        double maxit = 1e+4, double eps = 1e-04,
                        String warm = "lambda", bool strong = true) {
  // initialize
  int n = X_tilde.nrow();
  int p = X_tilde.ncol() + 1;
  int nlambda = lambda.size();

  // X
  NumericMatrix X_tilde2(n, p);
  for (int i = 0; i < n; i++) {
    X_tilde2(i, 0) = 1;
  }
  for (int j = 1; j < p; j++) {
    for (int i = 0; i < n; i++) {
      X_tilde2(i, j) = X_tilde(i, j - 1);
    }
  }

  // beta
  NumericMatrix beta(p, nlambda);
  NumericVector beta_tilde(p);
  NumericVector eta(n);

  // similarity matrix
  NumericMatrix R2(p, p);
  for (int j = 0; j < p; j++) {
    R2(0, j) = 0;
  }
  for (int j = 0; j < p; j++) {
    R2(0, j) = 0;
  }
  for (int j1 = 1; j1 < p; j1++) {
    for (int j2 = 1; j2 < p; j2++) {
      R2(j1, j2) = R(j1 - 1, j2 - 1);
    }
  }

  // candidate variables
  LogicalVector candid(p);
  candid(0) = true;
  NumericVector pi(n), w(n), s(n), r(n);

  // CDA for all lambda
  for (int k = 0; k < nlambda; k++) {
    // Rcout << "lambda " << k << std::endl;

    // beta-init or warm-start
    if (warm=="delta") {
      beta(_, k) = init_beta(_, k);
    } else if (warm=="lambda") {
      if (k != 0) {
        beta(_, k) = beta(_, k - 1);
      }
    }

    int iter = 0;
    while(iter < maxit) {
      while(iter < maxit) {
        iter++;

        beta_tilde = beta(_, k);

        for (int i = 0; i < n; i++) {
          eta(i) = 0;
          for (int j = 0; j < p; j++) {
            eta(i) = eta(i) + X_tilde2(i, j) * beta_tilde(j);
          }
          if (eta(i) > 10) {
            pi(i) = 1;
          } else if (eta(i) < -10) {
            pi(i) = 0;
          } else {
            pi(i) = 1 / (1 + exp(- eta(i)));
          }
          w(i) = pi(i) * (1 - pi(i));
          if (w(i) < 0.0001) {w(i) = 0.0001;}
          s(i) = y(i) - pi(i);
          r(i) = s(i) / w(i);
        }

        double maxChange = 0;

        for (int j = 0; j < p; j++) {
          if(candid(j)) {
            double xwr = 0;
            for (int i = 0; i < n; i++) {
              xwr = xwr + X_tilde2(i, j) * w(i) * r(i);
            }
            double xwx = 0;
            for (int i = 0; i < n; i++) {
              xwx = xwx + X_tilde2(i, j) * w(i) * X_tilde2(i, j);
            }

            double u, v, l1, l2;
            if (delta == 0) { // normal lasso
              u = xwr / n + (xwx / n) * beta_tilde(j);
              v = xwx / n;
              l1 = lambda(k);
              l2 = 0;
            } else { // exclusive lasso
              double rb = 0;
              for (int jj = 0; jj < p; jj++) {
                if ((jj != j) & (beta(jj, k) != 0)) {
                  rb = rb + R2(j, jj) * beta(jj, k);
                }
              }
              u = xwr / n + (xwx / n) * beta_tilde(j) - delta * lambda(k) * rb;
              v = xwx / n + delta * lambda(k) * R2(j, j);
              l1 = lambda(k);
              l2 = 0;
            }

            if (j == 0) {
              // beta(j, k) = updateLassoC(u, 0, l2, 1);
              // beta(j, k) = u / v;
              beta(j, k) = xwr / xwx + beta_tilde(j);
            } else {
              beta(j, k) = updateLassoC(u, l1, l2, v);
            }

            double shift = beta(j, k) - beta_tilde(j);
            if (shift != 0) {
              for (int i = 0; i < n; i++) {
                double si = shift * X_tilde2(i, j);
                r(i) = r(i) - si;
                eta(i) = eta(i) + si;
              }
            }

            if (maxChange < fabs(shift) * sqrt(v)) {
              maxChange = fabs(shift) * sqrt(v);
            }
          }
        }

        // check for convergence
        // double maxChange = max(abs(beta(_, k) - beta_tilde));
        beta_tilde = beta(_, k);
        if (maxChange < eps) {
          // Rcout << iter << std::endl;
          break;
        }
      }

      int violations = 0;
      for (int j = 0; j < p; j++) {
        if (!candid(j)) {
          double z = 0;
          for (int i = 0; i < n; i++) {
            z = z + X_tilde2(i, j) * s(i) / n;
          }
          if (fabs(z) > lambda(k)) {
            candid(j) = true;
            violations++;
          }
        }
      }
      if (violations == 0) {
        // Rcout << iter << std::endl;
        break;
      }
    }

    // break if maxit reached
    if (iter >= maxit) {
      return(beta);
    }
  }

  return beta;

}
