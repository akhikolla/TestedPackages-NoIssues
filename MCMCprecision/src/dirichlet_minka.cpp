#include <Rcpp.h>
using namespace Rcpp;

// Inverse digamma function (Minka, 2000)
// [[Rcpp::export]]
NumericVector inv_digamma(NumericVector y,
                          int iter = 5)
{
  NumericVector x = exp(y) + 0.5;
  double gamma = - R::psigamma(1., 0);
  for (int i = 0; i < y.size(); i++)
  {
    if (y[i] < -2.22)
      x[i] = -1/(y[i] + gamma);

    // # Newton iterations
    for(int k = 0; k < iter; k++)
      x[i] = x[i] - (R::digamma(x[i]) - y[i]) / R::psigamma(x[i], 1);
  }
  return x;
}

// Estimate dirichlet parameters with fixpoint-algorithm by Minka (2002)
// [[Rcpp::export]]
NumericVector dirichlet_fp(NumericVector alpha,
                           NumericVector logx_mean,
                           int maxit = 1e5,
                           double abstol = 1e-5)
{
  NumericVector alpha0 = alpha;
  int cnt = 0;
  double diff = 1.;

  while ((diff > abstol) && (cnt < maxit)) {
    alpha0 = alpha;
    // Fixpoint iteration:
    alpha = inv_digamma(R::digamma(sum(alpha0)) + logx_mean);

    cnt += 1;
    // alpha[alpha < 0] = min;
    diff = max(abs(alpha0 - alpha));
  }
  if (cnt == maxit)
    warning("Maximum number of iterations reached.");
  return alpha;
}

