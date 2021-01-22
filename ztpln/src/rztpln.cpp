#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>  // remove and remove_if
using namespace Rcpp;

// [[Rcpp::export]]
long do_rztpln(double mu, double sig) {
  double tau, lambda;
  tau = R::rnorm(0, 1);
  lambda = exp(mu + tau * sig);
  return R::qpois(R::runif(exp(-lambda), 1), lambda, 1, 0);
}

// [[Rcpp::export]]
long do_rpln(double mu, double sig) {
  double tau, lambda;
  tau = R::rnorm(0, 1);
  lambda = exp(mu + tau * sig);
  return R::qpois(R::runif(0, 1), lambda, 1, 0);
}

// [[Rcpp::export]]
Rcpp::IntegerVector do_vec_rztpln2(int n, double mu, double sig) {
  Rcpp::IntegerVector y(n);
  for (int i = 0; i < n; i++) {
    y(i)  = do_rztpln(mu, sig);
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::IntegerVector do_vec2_rztpln2(int n, Rcpp::NumericVector mu, 
    Rcpp::NumericVector sig) {
  Rcpp::IntegerVector y(n);
  for (int i = 0; i < n; i++) {
    y(i)  = do_rztpln(mu(i), sig(i));
  }
  return y;
}

/* repeat do_rpln and remove zero until it reached length n */

// [[Rcpp::export]]
Rcpp::IntegerVector do_vec_rztpln1(int n, double mu, double sig) {
  Rcpp::IntegerVector y(0);
  while (y.size() < n) {
    for (int i = y.size(); i < n; i++) {
     y.push_back(do_rpln(mu, sig));
    }
    y.erase(std::remove(y.begin(), y.end(), 0), y.end());
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::IntegerVector do_vec2_rztpln1(int n, Rcpp::NumericVector mu, 
    Rcpp::NumericVector sig) {
  Rcpp::IntegerVector y(0);
  while (y.size() < n) {
    for (int i = y.size(); i < n; i++) {
     y.push_back(do_rpln(mu(i), sig(i)));
    }
    y.erase(std::remove(y.begin(), y.end(), 0), y.end());
  }
  return y;
}
