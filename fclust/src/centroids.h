#ifndef CENTROIDS_H
#define CENTROIDS_H

#include <RcppArmadillo.h>

arma::mat centroids_FKM(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p, double m);

arma::mat centroids_FKM_ent(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p);

arma::mat centroids_FKM_pf(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p, double b);


#endif
