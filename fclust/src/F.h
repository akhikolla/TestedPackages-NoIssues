#ifndef F_H
#define F_H

#include <RcppArmadillo.h>

arma::cube F_gkb(arma::mat data, arma::mat U, arma::mat H, arma::mat F0, double m, double gam, unsigned int n, unsigned int k, int p, double mcn, arma::vec vp);

arma::cube F_gkb_ent(arma::mat data, arma::mat U, arma::mat H, arma::mat F0, double gam, unsigned int n, unsigned int k, int p, double mcn, arma::vec vp);

arma::cube F_gk(arma::mat data, arma::mat U, arma::mat H, double m,unsigned int n, unsigned int k, int p, arma::vec vp );

arma::cube F_gk_ent(arma::mat data, arma::mat U, arma::mat H,unsigned int n, unsigned int k, int p, arma::vec vp );
  

#endif
