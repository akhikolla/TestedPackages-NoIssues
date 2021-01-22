#ifndef __ProbitMain_h__
#define __ProbitMain_h__

#include "iBATCGH_RcppExports.h"

double MHR2(double phi,int* mR,double* mX,double* mY,double c,double d,double delta,double* mSum,arma::icolvec *Rout,double alpha0,double alpha1,int* mXI,arma::mat Hs, arma::vec countnotwo, int nR, int choiceg, int selectioncgh,int m,int mg, int g, int nOss);
double MHxi2(int* mXI,arma::mat A,arma::mat Picum,double* mDist,double disfix,double* mSum,int* mR,double* mY,double c,double d,double delta,arma::colvec mu,arma::colvec sigma,double* mX,arma::mat tran,arma::colvec Pi,arma::imat *xiout,arma::colvec *sumout,double den,double alpha0,double alpha1,arma::mat Hs,int choice,int sel,int m,int mg,int g,int nOss);
arma::vec sm(arma::mat distance,arma::imat xi,double disfix);
arma::colvec initial(arma::mat A);
arma::rowvec dir(arma::rowvec a);
double truncnorm(double mean,double sd,double lower,double upper);
double truncgamma(double a,double b,double lower);

#endif // __ProbitMain_h__

