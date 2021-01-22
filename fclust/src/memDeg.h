#ifndef MEMDEG_H
#define MEMDEG_H

#include <RcppArmadillo.h>

arma::mat memb_degree(arma::mat D, double m, unsigned int n, unsigned int k, unsigned int p);

arma::mat memb_degree_pf(arma::mat D, double b, unsigned int n, unsigned int k, unsigned int p);

arma::mat memb_degree_ent(arma::mat D, double ent, unsigned int n, unsigned int k, unsigned int p);

arma::mat memb_degree_medoid(arma::mat D, arma::uvec medoid, double m, unsigned int n, unsigned int k, unsigned int p);


#endif
