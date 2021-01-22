# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

using namespace Rcpp;

double ddc_ (double x, NumericVector phi) {
  int k = phi.length();
  arma::mat invB(k+1, k+1);
  arma::mat c(k+1, 1);
  arma::mat res(k, k);
  
  double normdens;
  
  // invB
  invB = invBMat(k);
  
  // Create c
  c = cMat(k, phi);
  
  res = arma::pow((invB * c).t() * expVec(x, k), 2);
  normdens = dnorm(NumericVector::create(x), 0, 1, true)[0];
  
  return exp(log(res[0]) + normdens);
}

//' Density function for univariate Davidian curves
//'
//' Returns the density for a vector of x.
//'
//' @param x vector of quantiles.
//' @param phi Davidian curve parameters.
//' A maximum of 10 parameters is allowed.
//' 
//' @examples
//' curve(ddc(x, 1.570789), -6, 6) # Approximately normal.
//' 
//' phi <- c(77.32, 78.51, 76.33, 77.16)
//' curve(ddc(x, phi), -6, 6) # A bimodal density.
//' integrate(ddc, phi = phi, lower = -Inf, upper = Inf) # Integrates to 1.
// [[Rcpp::export]]
NumericVector ddc (NumericVector x, NumericVector phi) {
  
  if (phi.length() > 10) {
    stop("length(phi) > 10 is not supported.");
  }

  NumericVector res(x.length());

  for (int i = 0; i < x.length(); i++) {
    if (Rcpp::traits::is_infinite<REALSXP>(x[i])) {
      res[i] = 0;
    } else {
      res[i] = ddc_(x[i], phi); 
    }
  }

  return res;
}

//' Random samples from univariate Davidian curves
//'
//' Returns n samples from a univariate Davidian curve.
//'
//' @param n Number of observations to be sampled.
//' @param phi Davidian curve parameters.
//' A maximum of 10 parameters is allowed.
//' 
//' @examples
//' # Sample from the standard normal Davidian curve:
//' hist(rdc(1000, 1.570789), xlim = c(-6, 6), ylim = c(0, 0.5), freq = FALSE, breaks = 20)
//' curve(dnorm(x), -6, 6, col = "blue", lwd = 1, add = TRUE)
//' curve(ddc(x, 1.570789), -6, 6, col = "red", lwd = 2, lty = 3, add = TRUE)
//'
//' # Sample from a bimodal density:
//' phi <- c(77.32, 78.51, 76.33, 77.16)
//' hist(rdc(1000, phi), xlim = c(-6, 6), ylim = c(0, 0.4), freq = FALSE, breaks = "fd")
//' curve(ddc(x, phi), -6, 6, col = "red", lwd = 2, lty = 3, add = TRUE)
//' 
// [[Rcpp::export]]
NumericVector rdc (int n, NumericVector phi) {
  
  if (phi.length() > 10) {
    stop("length(phi) > 10 is not supported.");
  }

  NumericVector out(n);
  NumericVector c(1, 8.55); // this may not be optimal but some local tests held up.
  int accepted = 0;
  
  NumericVector y(1);
  NumericVector u(1);
  NumericVector ratio(1);
  
  while (accepted < n) {
    y = runif(1, -10.0, 10.0); // proposal density is runif, domain is defined as [-10, +10].
    u = runif(1);
    
    ratio = ddc(y, phi) / (c * dunif(y, -10.0, +10.0)); //  ;
    
    if (is_true(all(u < ratio))) {
      out[accepted] = y[0];
      accepted++;
    }
  }
  
  return out;
}

