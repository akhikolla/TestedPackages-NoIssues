#ifndef __DCOV__
#define __DCOV__

#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;



int centering_from_data(const arma::mat& x,arma::mat &D,std::string type="U");
double dcov1v1(arma::vec x,arma::vec y,std::string type="U");
double dcor1v1(arma::vec x,arma::vec y,std::string type="U");
double dcor(const arma::mat &x,const arma::mat &y, std::string type="V");
double pdcov1v1v1(arma::vec x, arma::vec y, arma::vec z, std::string type);
double pdcor1v1v1(arma::vec x, arma::vec y, arma::vec z, std::string type);
double psum(arma::vec& bit, int i);
void update(arma::vec& bit, int n, int i, double val);


#endif
