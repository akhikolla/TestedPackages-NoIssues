#ifndef DISTANCE_H
#define DISTANCE_H

#include <RcppArmadillo.h>

arma::mat euclidean_distance(arma::mat data, arma::mat H, unsigned int n, unsigned int k, bool Square = false);

arma::mat euclidean_distance_gkb(arma::mat data, arma::mat H, arma::cube F, unsigned int n, unsigned int k, bool Square = false);

arma::mat euclidean_distance_gk(arma::mat data, arma::mat H, arma::cube F, arma::mat D_old, unsigned int n, unsigned int k, unsigned int p,bool Square = false);

arma::mat euclidean_distance_medoid(arma::mat data, arma::mat H, unsigned int n, unsigned int k, bool Square = false);


#endif
