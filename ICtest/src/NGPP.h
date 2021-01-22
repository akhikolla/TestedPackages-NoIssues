#include <RcppArmadillo.h>

RcppExport arma::mat symmetricPower_C(arma::mat x, double r);
RcppExport arma::vec computeObj_C(arma::mat x, arma::vec nl, arma::vec alpha);
RcppExport arma::vec computeTVec_C(arma::vec u, arma::mat x, arma::vec nl, arma::vec alpha);
