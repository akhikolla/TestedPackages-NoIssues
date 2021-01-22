#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List poly_divide(NumericVector Num0, NumericVector Den)
{
  int n = Num0.size(), m = Den.size();
  NumericVector Num = clone(Num0);
  if(n < m) {
    NumericVector Quo(1);
    return List::create(_["quotient"] = Quo,
                        _["remainder"] = Num);
  }
  int p = n - m;
  NumericVector Quo(p+1);
  for(int i = p; i >= 0; --i) {
    int ntop = (n-1)-(p-i);
    double f = Num(ntop)/Den(m-1);
    Quo(i) = f;
    for(int j = 0; j < m; j++) {
      Num(ntop - j) -= f * Den(m - 1 - j);
    }
    Num(ntop) = 0.0;
  }
  int rsize = m > 1 ? m - 1 : 1;
  NumericVector Rem(rsize);
  for(int j = 0; j < rsize; j++) {
    Rem(j) = Num(j);
  }
  return List::create(_["quotient"] = Quo,
                      _["remainder"] = Rem);
}

// [[Rcpp::export]]
NumericVector poly_product(NumericVector P, NumericVector Q)
{
  int np = P.size(), nq = Q.size();
  NumericVector PQ(np + nq - 1);
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < nq; j++) {
      PQ(i+j) += P(i)*Q(j);
    }
  }
  return PQ;
}

/*** R
poly_divide(poly_product(1:3, 5:1), 1:3)
poly_divide(poly_product(1:3, 5:1), 5:1)
poly_divide(poly_product(1:3, 5:1), 1:4)
poly_divide(poly_product(1:3, 5:1), 1:5)
*/
