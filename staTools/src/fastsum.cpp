#include <Rcpp.h>
using namespace Rcpp;
//' Fast Sum implemented in Cpp.
//'
//' Cpp function which speed up the computation of the Hurwitz zeta function.
//' @param i An integer.
//' @param xmin An integer.
//' @param alpha A real number greater than 1.
//' @export
//[[Rcpp::export]]

double fastsum(int i, int xmin, double alpha)
{
              double r = 0.0;
              double power = 0.0;
              for (double n=xmin; n<=i; n++) {
              power = pow(n, alpha);
              r = r + (1.0 / power );
              }
              return r;
}
