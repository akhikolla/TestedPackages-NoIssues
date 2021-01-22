#ifndef INCDTW_H
#define INCDTW_H

#include <RcppArmadillo.h>

double mymin(double x, double y);

Rcpp::NumericMatrix cpp_cm(const arma::mat& x,
                           const arma::mat& y,
                           std::string dist_method, int ws);

double dist_norm1(const arma::mat& x,
                  const arma::mat& y,
                  int i, int j, int ncol);

double dist_norm2(const arma::mat& x,
                  const arma::mat& y,
                  int i, int j, int ncol); 

double dist_norm2_square(const arma::mat& x,
                         const arma::mat& y,
                         int i, int j, int ncol);

typedef double (*funcPtr_dist)(const arma::mat& x,
                               const arma::mat& y,
                               int i, int j, int ncol);


typedef double (*funcPtr_dtw_mv)(const arma::mat& x, const arma::mat& y, 
                const Rcpp::List params);

double mystep_symmetric1 (const double gcm10, const double gcm11, const double gcm01, const double cm00);

double mystep_symmetric2 (const double gcm10, const double gcm11, const double gcm01, const double cm00);


double cpp_dtw2vec         (const arma::vec& x, const arma::vec& y, std::string step_pattern);
double cpp_dtw2vec_ws      (const arma::vec& x, const arma::vec& y, std::string step_pattern, int ws);
double cpp_dtw2vec_ea      (const arma::vec& x, const arma::vec& y, std::string step_pattern, double threshold);
double cpp_dtw2vec_ws_ea   (const arma::vec& x, const arma::vec& y, std::string step_pattern, int ws, double threshold);

double cpp_dtw2vec_v32     (const arma::vec& x, const arma::vec& y);
double cpp_dtw2vec_dummy   (const arma::vec& x, const arma::vec& y, std::string dummy);
   
double cpp_dtw2vec_mv      (const arma::mat& x, const arma::mat& y, const Rcpp::List params);
double cpp_dtw2vec_mv_ws_ea(const arma::mat& x, const arma::mat& y, std::string step_pattern, std::string dist_method, int ws, double threshold);

double multp_dtw2vec_ws_ea   (const arma::vec& x, const arma::vec& y, std::string step_pattern, int ws, double threshold);
double multp_dtw2vec_mv_ws_ea(const arma::mat& x, const arma::mat& y, std::string step_pattern, std::string dist_method, int ws, double threshold);



#endif
