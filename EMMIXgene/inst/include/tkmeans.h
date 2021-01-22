#ifndef PACKAGENAME_TKMEANS_H
#define PACKAGENAME_TKMEANS_H

arma::mat tkmeans(arma::mat& M, int k , double alpha, arma::vec weights,  int nstart = 1, int iter = 10, double tol = 0.0001, bool verbose = false);

#endif
