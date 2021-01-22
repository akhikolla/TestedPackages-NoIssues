#ifndef FUNC_H
#define FUNC_H

#include <RcppArmadillo.h>


arma::rowvec replace(arma::rowvec x, double val, double change);

arma::mat InvCheck(arma::mat A);

bool Match(int i, arma::uvec B);

void distCheck(Rcpp::NumericMatrix D, unsigned int n, unsigned int p);

double silhouette_internal(arma::mat X, arma::mat U, unsigned int p, unsigned int k, unsigned int n,bool distance = false);

Rcpp::List silhouette(arma::mat X, arma::mat U, unsigned int p, unsigned int k, unsigned int n, bool distance = false);

double silhouetteFuzzy(arma::mat X, arma::mat U, double alpha,unsigned int p, unsigned int k, unsigned int n, bool distance = false);

double partCoef(arma::mat U, unsigned int n);

double partEntropy(arma::mat U, double b, unsigned int n);

double partCoef_mod(arma::mat U, unsigned int n, unsigned int k);

double xie_beni(arma::mat X,arma::mat U, arma::mat H, double m, unsigned int n, unsigned int k);

double indices(std::string type, arma::mat X, arma::mat U, arma::mat H, double m, unsigned int n, unsigned int k, unsigned int p, double b, double alpha, bool distance = false);

arma::mat unifInit(int n, int d);
  

#endif


