#include <RcppArmadillo.h>
// using namespace Rcpp;
using namespace std;

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double mymin(double x, double y)
{
   // z > nan for z != nan is required by C the standard
   int xnan = isnan(x), ynan = isnan(y);
   if(xnan || ynan) {
      if(xnan && !ynan) return y;
      if(!xnan && ynan) return x;
      return x;
   }
   return std::min(x,y);
}



// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist_norm1(const arma::mat& x,
                  const arma::mat& y,
                  int i, int j, int ncol) {
   double ret = 0;
   for(int k =0; k < ncol; k++){
      ret += abs( x(i, k) - y(j, k) );
   }
   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist_norm2(const arma::mat& x,
                  const arma::mat& y,
                  int i, int j, int ncol) {
   double ret = 0;
   double tmp = 0;
   for(int k =0; k < ncol; k++){
      tmp = x(i, k) - y(j, k);
      ret += tmp * tmp;
   }
   ret = sqrt(ret);
   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double dist_norm2_square(const arma::mat& x,
                         const arma::mat& y,
                         int i, int j, int ncol) {
   double ret = 0;
   double tmp = 0;
   for(int k =0; k < ncol; k++){
      tmp = x(i, k) - y(j, k);
      ret += tmp * tmp;
   }
   ret = ret * 1;
   return (ret);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*funcPtr_dist)(const arma::mat& x,
                               const arma::mat& y,
                               int i, int j, int ncol);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double mystep_symmetric1(const double gcm10, //vertical
                         const double gcm11, //diagonal
                         const double gcm01, //horizontal
                         const double cm00){
   return(cm00 + mymin(gcm10, mymin(gcm11, gcm01)));
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double mystep_symmetric2(const double gcm10, //vertical
                         const double gcm11, //diagonal
                         const double gcm01, //horizontal
                         const double cm00){
   return(cm00 + mymin(gcm10, mymin(cm00 + gcm11, gcm01)));
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*funcPtr_step)(const double gcm10, 
                               const double gcm11, 
                               const double gcm01, 
                               const double cm00);


