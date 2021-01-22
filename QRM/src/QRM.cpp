#include <Rcpp.h>
using namespace Rcpp;

double contribH(double s, double y, double gammaval, double deltaval)
{
  double result;
  result = (1+deltaval*y)*exp(-gammaval*s);
  return result;
}

double contribE(double s, double y, double gammaval, double rhoval, double deltaval)
{
  double result;
  result = (1+deltaval*y)/pow(1+s/gammaval,(1.0+rhoval));
  return result;
}

double ef(double x, double lambda, double chi, double psi) {
  double result;
  result = pow(x, lambda - 1.0) * exp(-0.5 * (psi * x + chi / x));
  return result;
}

// [[Rcpp::export]]
NumericVector rfrank(int n, double theta) {

  int i, k;
  double p = 0, alpha;
  NumericVector ans(n), U(n);

  alpha = 1 - exp(-theta);
  U = runif(n, 0, 1);

  for(i = 0; i < n; i++){
    k = 1;
    p = alpha / theta;
    while(U[i] > p){
      k++;
      p = p + pow(alpha, k) / (k * theta);
    }
    ans[i] = k;
  }
  return ans;
}

// [[Rcpp::export]]
NumericVector rgig(int n, double r, double s, double p, double k1, double k2, double lambda, double chi, double psi, double s1, double s2) {

  int i = 0;
  double U, Ustar, level, x;
  NumericVector ans(n);

  while (i < n){
    U = R::runif(0, 1);
    Ustar = R::runif(0, 1);

    if (U <= r){
      x = log(1 + s * U / k1) / s;  
      level = log(ef(x, lambda, chi, psi + 2 * s) / s1);
      if (log(Ustar) <= level){
	ans[i] = x;
	i++;
      }
    } else {
      x = -log(p * (1 - U) / k2) / p;
      level = log(ef(x, lambda, chi, psi - 2 * p) / s2);
      if (log(Ustar) <= level){
	ans[i] = x;
	i++;
      }
    }
  }

  return ans;
}

